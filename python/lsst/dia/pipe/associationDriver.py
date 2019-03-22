import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet

from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.base import ArgumentParser, TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError
from lsst.pipe.drivers.utils import getDataRef, TractDataIdContainer
from lsst.pipe.tasks.coaddBase import getSkyInfo

from .simple_association import SimpleAssociationTask

__all__ = ("AssociationDriverConfig", "AssociationDriverTask")


class AssociationDriverConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")
    maxFootprintArea = Field(dtype=int, default=5000, doc="Maximum area of footprints")
    defaultFootprintRadius = Field(
        dtype=int, default=40,
        doc="Use this radius to define the footprint if too large"
    )
    associator = ConfigurableField(
        target=SimpleAssociationTask,
        doc="Task used to associate DiaSources with DiaObjects.",
    )
    ccdKey = Field(dtype=str, default='detector', doc="Name of ccd to give to butler")
    maxLoadCatalogs = Field(
        dtype=int, default=None,
        doc="Use this radius to define the footprint if too large"
    )


class AssociationDriverTaskRunner(TaskRunner):
    """TaskRunner for running AssociationDriver Task.

    This is similar to the lsst.pipe.drivers.multiBandDriver runner
    except with no reuse outputs option.
    """

    def __init__(self, TaskClass, parsedCmd, doReturnResults=False):
        TaskRunner.__init__(self, TaskClass, parsedCmd, doReturnResults)

    def makeTask(self, parsedCmd=None, args=None):
        """A variant of the base version that passes a butler argument to the task's constructor
        parsedCmd or args must be specified.
        """
        if parsedCmd is not None:
            butler = parsedCmd.butler
        elif args is not None:
            dataRefList, kwargs = args
            butler = dataRefList[0].butlerSubset.butler
        else:
            raise RuntimeError("parsedCmd or args must be specified")
        return self.TaskClass(config=self.config, log=self.log, butler=butler)


def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)


def _makeGetSchemaCatalogs(datasetSuffix):
    """Construct a getSchemaCatalogs instance method
    datasetSuffix:  Suffix of dataset name, e.g., "src" for "deepCoadd_src"
    """

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        src = afwTable.SourceCatalog(self.schema)
        if hasattr(self, "algMetadata"):
            src.getTable().setMetadata(self.algMetadata)
        return {self.config.coaddName + "Diff_" + datasetSuffix: src}
    return getSchemaCatalogs


class AssociationDriverTask(BatchPoolTask):
    """Multi-node driver for multiband processing"""
    ConfigClass = AssociationDriverConfig
    _DefaultName = "associationDriver"
    RunnerClass = AssociationDriverTaskRunner
    getSchemaCatalogs = _makeGetSchemaCatalogs("diaObject")

    def __init__(self, butler=None, **kwargs):
        """!
        @param[in] butler: the butler can be used to retrieve schema or passed to the refObjLoader
            constructor in case it is needed.
        """
        BatchPoolTask.__init__(self, **kwargs)
        self.butler = butler
        self.src_schema = self.butler.get(f"{self.config.coaddName}Diff_diaSrc_schema").schema
        self.makeSubtask("associator", src_schema=self.src_schema)
        self.schema = afwTable.SourceTable.makeMinimalSchema()

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                                                   parentTask=self._parentTask, log=self.log,
                                                   butler=self.butler))

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name=cls._DefaultName, *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd",
                               help="data ID, e.g. --id tract=12345 patch=1,2 filter=g^r^i",
                               ContainerClass=TractDataIdContainer)
        return parser

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numCpus):
        """!Return walltime request for batch job

        @param time: Requested time per iteration
        @param parsedCmd: Results of argument parsing
        @param numCores: Number of cores
        """
        numTargets = 0
        for refList in parsedCmd.id.refList:
            numTargets += len(refList)
        return time*numTargets/float(numCpus)

    @abortOnError
    def runDataRef(self, patchRefList):
        """!Run association processing on coadds

        Only the master node runs this method.

        @param patchRefList:  Data references to run measurement
        """
        for patchRef in patchRefList:
            if patchRef:
                butler = patchRef.getButler()
                break
        else:
            raise RuntimeError("No valid patches")
        pool = Pool("all")
        pool.cacheClear()
        pool.storeSet(butler=butler)

        # Group all filters by patch
        patches = {}
        tract = None
        for patchRef in patchRefList:
            dataId = patchRef.dataId
            if tract is None:
                tract = dataId["tract"]
            else:
                assert tract == dataId["tract"]

            patch = dataId["patch"]
            if patch not in patches:
                patches[patch] = []
            patches[patch].append(dataId)

        pool.map(self.runAssociation, patches.values())

    def runAssociation(self, cache, dataIdList):
        """! Run association on a patch
        For all of the visits that overlap this patch in the band create a DIAObject
        catalog.  Only the objects in the non-overlaping area of the tract and patch
        are included.
        """

        dataRefList = [getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd_calexp") for
                       dataId in dataIdList]

        idFactory = None

        for dataRef in dataRefList:

            try:
                calexp = dataRef.get(f"{self.config.coaddName}Coadd_calexp")
            except Exception:
                self.log.info('Cannot read data for %s' % (dataRef.dataId))
                continue

            coaddWcs = calexp.getWcs()
            visitCatalog = calexp.getInfo().getCoaddInputs().ccds
            band = dataRef.dataId['filter']
            skyInfo = getSkyInfo(coaddName=self.config.coaddName, patchRef=dataRef)

            innerPatchBox = afwGeom.Box2D(skyInfo.patchInfo.getInnerBBox())
            self.log.info('Total number of images from %s %s %s to read %d' %
                          (dataRef.dataId['tract'], dataRef.dataId['patch'], band, len(visitCatalog)))

            numCatalogsLoaded = 0
            for visitRec in visitCatalog:

                if numCatalogsLoaded >= self.config.maxLoadCatalogs and self.config.maxLoadCatalogs:
                    continue
                visit = visitRec.get('visit')
                ccd = visitRec.get('ccd')
                dataId = {"visit": visit, self.config.ccdKey: ccd}
                try:
                    exp = cache.butler.get(f"{self.config.coaddName}Diff_differenceExp", dataId)
                    src = cache.butler.get(f"{self.config.coaddName}Diff_diaSrc", dataId)
                    srcWcs = exp.getWcs()
                except Exception as e:
                    self.log.debug('Cannot read data for %d %d. skipping %s' % (visit, ccd, e))
                    continue
                numCatalogsLoaded += 1
                if idFactory is None:
                    expBits = dataRef.get("deepMergedCoaddId_bits")
                    expId = int(dataRef.get("deepMergedCoaddId"))
                    idFactory = afwTable.IdFactory.makeSource(expId, 64 - expBits)
                    self.associator.initialize(idFactory)

                mask = np.array(
                    [innerPatchBox.contains(coaddWcs.skyToPixel(srcWcs.pixelToSky(a.getCentroid())))
                     for a in src],
                    dtype=bool
                )

                src = src[mask]
                if len(src) == 0:
                    continue

                self.log.info('Reading difference image %d %d, %s with %d possible sources' %
                              (visit, ccd, dataRef.dataId['filter'], len(src)))

                footprints = []
                region = calexp.getBBox(afwImage.PARENT)
                for ii, rec in enumerate(src):
                    # transformations on large footprints can take a long time
                    # We truncate the footprint since we will rarely be interested
                    # in such large footprints
                    if rec.getFootprint().getArea() > self.config.maxFootprintArea:
                        spans = afwGeom.SpanSet.fromShape(self.config.defaultFootprintRadius,
                                                          afwGeom.Stencil.CIRCLE,
                                                          afwGeom.Point2I(rec.getCentroid()))
                        foot = afwDet.Footprint(spans)
                        foot.addPeak(int(rec.getX()), int(rec.getY()), 1)
                    else:
                        foot = rec.getFootprint()
                    footprints.append(foot.transform(srcWcs, coaddWcs, region))

                self.associator.addCatalog(src, band, visit, ccd, footprints)

        result = self.associator.finalize(idFactory)

        if len(dataRefList) > 0 and result is not None:
            dataRefList[0].put(result, self.config.coaddName + 'Diff_diaObject')
            self.log.info('Total objects found %d' % len(result))

    def writeMetadata(self, dataRef):
        """We don't collect any metadata, so skip"""
        pass

    def _getConfigName(self):
        """!Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return "associationDriver_config"

    def _getMetadataName(self):
        """!Return the name of the metadata dataset.  Forced metadata to be saved
        """
        return "associationDriver_metadata"

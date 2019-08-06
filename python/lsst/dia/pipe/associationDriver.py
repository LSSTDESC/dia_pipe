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
from lsst.pipe.tasks.coaddBase import getSkyInfo, WcsSelectImagesTask
from lsst.pipe.tasks.coaddBase import SelectDataIdContainer
from lsst.pipe.base import Struct

from .simple_association import SimpleAssociationTask

__all__ = ("AssociationDriverConfig", "AssociationDriverTask")


class AssociationDriverConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")
    maxFootprintArea = Field(dtype=int, default=5000, doc="Maximum area of footprints")
    defaultFootprintRadius = Field(
        dtype=int, default=40,
        doc="Use this radius to define the footprint if too large")
    associator = ConfigurableField(
        target=SimpleAssociationTask,
        doc="Task used to associate DiaSources with DiaObjects.",
    )
    ccdKey = Field(dtype=str, default='detector', doc="Name of ccd to give to butler")
    select = ConfigurableField(
        doc="Image selection subtask.",
        target=WcsSelectImagesTask,
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

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return TaskRunner.getTargetList(parsedCmd, selectDataList=parsedCmd.selectId.dataList, **kwargs)


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
        self.makeSubtask("associator")
        self.makeSubtask("select")
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
        parser.add_id_argument("--selectId", "deepDiff_differenceExp",
                               help="data ID, e.g. --selectId visit=6789 ccd=0..9",
                               ContainerClass=SelectDataIdContainer)
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
    def runDataRef(self, patchRefList, selectDataList=[]):
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

        pool.map(self.runAssociation, patches.values(), selectDataList)

    def catalogGenerator(self, cache, dataRefList, selectDataList=[]):
        """! Get a generator of difference images from the visit catalof for each datRef
        """
        for dataRef in dataRefList:
            try:
                calexp = dataRef.get(f"{self.config.coaddName}Coadd_calexp")
            except Exception:
                self.log.info('Cannot coadd read data for %s' % (dataRef.dataId))
                continue

            visitCatalog = calexp.getInfo().getCoaddInputs().ccds
            for visitRec in visitCatalog:

                visit = visitRec.get('visit')
                ccd = visitRec.get('ccd')
                dataId = {"visit": visit, self.config.ccdKey: ccd}
                try:
                    exp = cache.butler.get(f"{self.config.coaddName}Diff_differenceExp", dataId)
                    src = cache.butler.get(f"{self.config.coaddName}Diff_diaSrc", dataId)
                except Exception as e:
                    self.log.debug('Cannot read difference data for %d %d. skipping %s' % (visit, ccd, e))
                    continue
                data = Struct(visit=visit, ccd=ccd, exp=exp, src=src,
                              filter=dataRef.dataId['filter'], calib=exp.getPhotoCalib())
                yield data

    def idListGenerator(self, cache, dataRefList, selectDataList=[]):
        """! Get a generator of difference images from the selectDataList
        """
        skyInfo = getSkyInfo(coaddName=self.config.coaddName, patchRef=dataRefList[0])
        diffExpRefList = self.selectExposures(dataRefList[0], skyInfo, selectDataList=selectDataList)
        for dataRef in diffExpRefList:
            try:
                exp = dataRef.get(f"{self.config.coaddName}Diff_differenceExp")
                src = dataRef.get(f"{self.config.coaddName}Diff_diaSrc")
            except Exception as e:
                self.log.debug('Cannot read data for %d %d. skipping %s' % (dataRef.dataId['visit'],
                                                                            dataRef.dataId['ccd'], e))
                continue
            data = Struct(visit=dataRef.dataId['visit'],
                          ccd=dataRef.dataId[self.config.ccdKey],
                          filter=dataRef.dataId['filter'], exp=exp, src=src, calib=exp.getPhotoCalib())
            yield data

    def runAssociation(self, cache, dataIdList, selectDataList):
        """! Run association on a patch
        For all of the visits that overlap this patch in the band create a DIAObject
        catalog.  Only the objects in the non-overlaping area of the tract and patch
        are included.
        """
        dataRefList = [getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd_calexp") for
                       dataId in dataIdList]

        # We need the WCS for the patch, so we can use the first entry in the dataIdList
        dataRef = dataRefList[0]
        skyInfo = getSkyInfo(coaddName=self.config.coaddName, patchRef=dataRef)
        try:
            calexp = dataRef.get(f"{self.config.coaddName}Coadd_calexp")
        except Exception:
            self.log.info('Cannot read coadd data for %s' % (dataRef.dataId))
            return

        coaddWcs = calexp.getWcs()
        innerPatchBox = afwGeom.Box2D(skyInfo.patchInfo.getInnerBBox())

        expBits = dataRef.get("deepMergedCoaddId_bits")
        expId = int(dataRef.get("deepMergedCoaddId"))
        idFactory = afwTable.IdFactory.makeSource(expId, 64 - expBits)

        if len(selectDataList) == 0:
            differenceImages = self.catalogGenerator
        else:
            differenceImages = self.idListGenerator

        initializeSelector = False
        for diffIm in differenceImages(cache, dataRefList, selectDataList):

            if initializeSelector is False:
                self.associator.initialize(diffIm.src.schema, idFactory)
                initializeSelector = True

            if len(diffIm.src) == 0:
                continue

            srcWcs = diffIm.exp.getWcs()
            isInside = np.array(
                [innerPatchBox.contains(coaddWcs.skyToPixel(srcWcs.pixelToSky(a.getCentroid())))
                 for a in diffIm.src],
                dtype=bool
            )

            isGood = np.array(
                [rec.getFootprint().contains(afwGeom.Point2I(rec.getCentroid()))
                 for rec in diffIm.src],
            )

            mask = (isInside) & (isGood)

            src = diffIm.src[mask]
            if len(src) == 0:
                continue

            self.log.info('Reading difference image %d %d, %s with %d possible sources' %
                          (diffIm.visit, diffIm.ccd, diffIm.filter, len(src)))

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

            self.associator.addCatalog(src, diffIm.filter, diffIm.visit, diffIm.ccd, diffIm.calib, footprints)

        result = self.associator.finalize(idFactory)

        if len(dataRefList) > 0 and result is not None:
            dataRefList[0].put(result, self.config.coaddName + 'Diff_diaObject')
            self.log.info('Total objects found %d' % len(result))

            idCatalog = self.associator.getObjectIds()
            dataRefList[0].put(idCatalog, self.config.coaddName + 'Diff_diaObjectId')

    def selectExposures(self, patchRef, skyInfo=None, selectDataList=[]):
        """!
        @brief Select exposures to associate
        Get the corners of the bbox supplied in skyInfo using @ref afwGeom.Box2D and convert the pixel
        positions of the bbox corners to sky coordinates using @ref skyInfo.wcs.pixelToSky. Use the
        @ref WcsSelectImagesTask_ "WcsSelectImagesTask" to select exposures that lie inside the patch
        indicated by the dataRef.
        @param[in] patchRef  data reference for sky map patch. Must include keys "tract", "patch",
                             plus the camera-specific filter key (e.g. "filter" or "band")
        @param[in] skyInfo   geometry for the patch; output from getSkyInfo
        @return    a list of science exposures to coadd, as butler data references
        """
        if skyInfo is None:
            skyInfo = self.getSkyInfo(patchRef)
        cornerPosList = afwGeom.Box2D(skyInfo.bbox).getCorners()
        coordList = [skyInfo.wcs.pixelToSky(pos) for pos in cornerPosList]
        return self.select.runDataRef(patchRef, coordList, selectDataList=selectDataList).dataRefList

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

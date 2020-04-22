import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.meas.base.forcedPhotCoadd import ForcedPhotCoaddTask, ForcedPhotCoaddConfig, ForcedPhotCoaddRunner

from .forcedPhotDia import DiaReferencesTask

__all__ = ("ForcedPhotCoaddDiaConfig", "ForcedPhotCoaddDiaTask")

class ForcedPhotCoaddDiaConfig(ForcedPhotCoaddConfig):
    skipMissing = pexConfig.Field(dtype=bool, default=True,
                                  doc="Skip getting references if they do not exist?")
    coaddName = pexConfig.Field(dtype=str, default='deep',
                                doc="Name of coadd")

    def setDefaults(self):
        ForcedPhotCoaddTask.ConfigClass.setDefaults(self)
        self.footprintDatasetName = 'diaObject'
        self.references.retarget(DiaReferencesTask)
        self.measurement.copyColumns = {"id": "dia_object_id", "coord_ra": "coord_ra",
                                        "coord_dec": "coord_dec"}
        self.measurement.plugins.names = ['base_DiaTransformedCentroid', 'base_PsfFlux']
        self.measurement.slots.centroid = 'base_DiaTransformedCentroid'
        self.measurement.slots.shape = None
        self.measurement.slots.apFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.psfFlux = None
        self.measurement.slots.calibFlux = None
        self.measurement.doReplaceWithNoise = False


class ForcedPhotCoaddDiaTask(ForcedPhotCoaddTask):
    """!
    A command-line driver for performing forced measurement on Coadd images from DIAObject catalogs.

    This task is a subclass of ForcedPhotCoaddTask, although it does not use most of the
    functionality defined there.  In this task we delegate the process of looking up
    the reference catalog inside the run function.  We allow for the references to come
    from multiple tracts and perform forced phometry on each tract and then combine them.
    The refernces are all taken from the inner patch boundary
    """

    ConfigClass = ForcedPhotCoaddDiaConfig
    RunnerClass = ForcedPhotCoaddRunner
    _DefaultName = "forcedPhotCoaddDia"
    dataPrefix = "deepCoadd_"

    def writeOutput(self, dataRef, sources):
        """!Write source table
        @param dataRef  Data reference from butler
        @param sources  SourceCatalog to save
        """
        dataRef.put(sources, "deepDiff_forced_diaObject",
                    flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

    def _getConfigName(self):
        """!Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return "forcedPhotCoaddDia_config"

    def _getMetadataName(self):
        """!Return the name of the metadata dataset.  Forced metadata to be saved
        """
        return None

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID with coadd keys"
                               "e.g. --id tract=1234 patch=1,2 filter=r^i")
        parser.add_argument("--psfCache", type=int, default=100, help="Size of CoaddPsf cache")
        return parser

    def runDataRef(self, dataRef, psfCache=None):
        """!Measure a single exposure using forced detection for a reference catalog.
        @param[in]  dataRef   An lsst.daf.persistence.ButlerDataRef
        @param[in]  psfCache  Size of PSF cache, or None. The size of the PSF cache can have
                              a significant effect upon the runtime for complicated PSF models.
        """
        refWcs = self.references.getWcs(dataRef)
        exposure = self.getExposure(dataRef)

        if psfCache is not None:
            exposure.getPsf().setCacheCapacity(psfCache)

        refCat = self.getReferences(dataRef, exposure)

        measCat = self.measurement.generateMeasCat(exposure, refCat, refWcs,
                                                   idFactory=self.makeIdFactory(dataRef))
        self.log.info("Performing forced measurement on %s" % (dataRef.dataId,))
        self.attachFootprints(measCat, refCat, exposure, refWcs, dataRef)

        exposureId = self.getExposureId(dataRef)
        self.measurement.run(measCat, exposure, refCat, refWcs, exposureId=exposureId)

        if self.config.doApCorr:
            self.applyApCorr.run(
                catalog=measCat,
                apCorrMap=exposure.getInfo().getApCorrMap()
            )
        self.catalogCalculation.run(measCat)

        self.writeOutput(dataRef, measCat)

    def getReferences(self, dataRef, exposure):
        """Return an iterable of reference sources which overlap the exposure

        Copied from forcedPhotCoaddTask, modief to call fetchInPactches from this class

        @param dataRef       Data reference from butler corresponding to the image to be measured;
                             should have tract, patch, and filter keys.
        @param exposure      lsst.afw.image.Exposure to be measured (not used by this implementation)
        All work is delegated to the references subtask; see CoaddSrcReferencesTask for information
        about the default behavior.
        """
        skyMap = dataRef.get(self.dataPrefix + "skyMap", immediate=True)
        tractInfo = skyMap[dataRef.dataId["tract"]]
        patch = tuple(int(v) for v in dataRef.dataId["patch"].split(","))
        patchInfo = tractInfo.getPatchInfo(patch)
        references = afwTable.SourceCatalog(self.references.schema)
        references.extend(self.fetchInPatches(dataRef, patchList=[patchInfo]))
        return references

    def fetchInPatches(self, dataRef, patchList):
        """!
        Copied from CoaddSrcReferencesTask and modified to allow loading deepDiff_diaObjects.

        The given dataRef must include the tract in its dataId.
        """
        dataset = "deepDiff_diaObject"
        tract = dataRef.dataId["tract"]
        butler = dataRef.butlerSubset.butler
        for patch in patchList:
            dataId = {'tract': tract, 'patch': "%d,%d" % patch.getIndex()}

            if not butler.datasetExists(dataset, dataId):
                if self.config.skipMissing:
                    continue
                raise pipeBase.TaskError("Reference %s doesn't exist" % (dataId,))
            self.log.info("Getting references in %s" % (dataId,))
            catalog = butler.get(dataset, dataId, immediate=True)

            for source in catalog:
                yield source

    def attachFootprints(self, sources, refCat, exposure, refWcs, dataRef):
        """Just use the footprints from the diaObject catalog
        """
        for refRecord, srcRecord in zip(refCat, sources):
            srcRecord.setFootprint(refRecord.getFootprint())



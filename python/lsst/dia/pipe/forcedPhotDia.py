import numpy as np

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwGeom
import lsst.geom as geom
import lsst.sphgeom as sphgeom

from lsst.meas.base.references import MultiBandReferencesConfig, BaseReferencesTask, CoaddSrcReferencesTask
from lsst.meas.base.forcedPhotCcd import ForcedPhotCcdTask, ForcedPhotCcdConfig
from lsst.meas.base.forcedPhotCoadd import ForcedPhotCoaddTask, ForcedPhotCoaddConfig, ForcedPhotCoaddRunner
from lsst.meas.base.forcedMeasurement import ForcedPluginConfig, ForcedPlugin
from lsst.meas.base.pluginRegistry import register


__all__ = ("ForcedPhotCcdDiaConfig", "ForcedPhotCcdDiaTask",
           "ForcedPhotCoaddDiaConfig", "ForcedPhotCoaddDiaTask")


class DiaReferencesConfig(CoaddSrcReferencesTask.ConfigClass):
    def validate(self):
        if self.filter is not None:
            raise pexConfig.FieldValidationError(
                field=MultiBandReferencesConfig.filter,
                config=self,
                msg="Filter should not be set for the multiband processing scheme")
        # Delegate to ultimate base class, because the direct one has a check we don't want.
        BaseReferencesTask.ConfigClass.validate(self)


class DiaReferencesTask(CoaddSrcReferencesTask):
    """!
    Dummy task so that the reference schema will match

    We do not get the reference tracts/patches using this class
    """
    ConfigClass = DiaReferencesConfig
    datasetSuffix = "diaObject"

    def __init__(self, butler=None, schema=None, **kwargs):
        """! Initialize the task.
        We cannot use the CoaddSrcReferencesTask init because of the hard coded "Coadd"

        Additional keyword arguments (forwarded to BaseReferencesTask.__init__):
         - schema: the schema of the detection catalogs used as input to this one
         - butler: a butler used to read the input schema from disk, if schema is None
        The task will set its own self.schema attribute to the schema of the output merged catalog.
        """
        BaseReferencesTask.__init__(self, butler=butler, schema=schema, **kwargs)
        if schema is None:
            assert butler is not None, "No butler nor schema provided"
            schema = butler.get("{}Diff_{}_schema".format(self.config.coaddName, self.datasetSuffix),
                                immediate=True).getSchema()
        self.schema = schema


class ForcedDiaTransformedCentroidConfig(ForcedPluginConfig):
    pass


@register("base_DiaTransformedCentroid")
class ForcedDiaTransformedCentroidPlugin(ForcedPlugin):
    """A centroid pseudo-algorithm for forced measurement that simply transforms sky coordinate
    from the reference catalog to the measurement coordinate system.  This does assumes that
    there is no centroid column in the refernce catalog.
    """

    ConfigClass = ForcedDiaTransformedCentroidConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        # Allocate x and y fields, join these into a single FunctorKey for ease-of-use.
        xKey = schema.addField(name + "_x", type="D", doc="transformed reference centroid column",
                               units="pixel")
        yKey = schema.addField(name + "_y", type="D", doc="transformed reference centroid row",
                               units="pixel")
        self.centroidKey = afwTable.Point2DKey(xKey, yKey)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        targetWcs = exposure.getWcs()
        targetPos = targetWcs.skyToPixel(refRecord.getCoord())
        measRecord.set(self.centroidKey, targetPos)


class ForcedPhotCcdDiaConfig(ForcedPhotCcdConfig):
    skipMissing = pexConfig.Field(dtype=bool, default=True,
                                  doc="Skip getting references if they do not exist?")
    coaddName = pexConfig.Field(dtype=str, default='deep',
                                doc="Name of coadd")

    def setDefaults(self):
        ForcedPhotCcdTask.ConfigClass.setDefaults(self)
        self.references.retarget(DiaReferencesTask)
        self.measurement.copyColumns = {"id": "dia_object_id", "coord_ra": "coord_ra",
                                        "coord_dec": "coord_dec"}
        self.measurement.plugins.names = ['base_SdssShape', 'base_DiaTransformedCentroid',
                                          'base_PsfFlux', 'base_LocalBackground',
                                          'base_PixelFlags']
        self.measurement.slots.centroid = 'base_DiaTransformedCentroid'
        self.measurement.slots.shape = 'base_SdssShape'
        self.measurement.slots.apFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.psfFlux = None
        self.measurement.slots.calibFlux = None


class ForcedPhotCcdDiaTask(ForcedPhotCcdTask):
    """!
    A command-line driver for performing forced measurement on CCD images from DIAObject catalogs.

    This task is a subclass of ForcedPhotCcdTask, although it does not use most of the
    functionality defined there.  In this task we delegate the process of looking up
    the reference catalog inside the run function.  We allow for the references to come
    from multiple tracts and perform forced phometry on each tract and then combine them.
    The refernces are all taken from the inner patch boundary
    """

    ConfigClass = ForcedPhotCcdDiaConfig
    RunnerClass = pipeBase.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCcdDia"
    dataPrefix = ""

    def writeOutput(self, dataRef, sources):
        """!Write source table
        @param dataRef  Data reference from butler
        @param sources  SourceCatalog to save
        """
        dataRef.put(sources, self.dataPrefix + "deepDiff_forced_diaSrc",
                    flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

    def _getConfigName(self):
        """!Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return self.dataPrefix + "forcedPhotCcdDia_config"

    def _getMetadataName(self):
        """!Return the name of the metadata dataset.  Forced metadata to be saved
        """
        return None

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepDiff_differenceExp", help="data ID with raw CCD keys"
                               "e.g. --id visit=12345 ccd")
        return parser

    def runDataRef(self, dataRef, psfCache=None):
        """!Measure a single exposure using forced detection for a reference catalog.
        @param[in]  dataRef   An lsst.daf.persistence.ButlerDataRef. It is passed to the
                              references subtask to obtain the reference WCS, the getExposure()
                              method (implemented by derived classes) to read the measurement
                              image, and load the reference catalog.    Sources are
                              generated with generateMeasCat() in the measurement subtask.  These
                              are passed to measurement's run method which fills the source
                              catalog with the forced measurement results.  The sources are then
                              passed to the writeOutputs() method (implemented by derived classes)
                              which writes the outputs.  See derived class documentation for which
                              datasets and data ID keys are used.
        @param[in]  psfCache  Size of PSF cache, or None. The size of the PSF cache can have
                              a significant effect upon the runtime for complicated PSF models.
        """
        exposure = self.getExposure(dataRef)
        tractRefCat = self.getReferences(dataRef.getButler(), exposure)

        if psfCache is not None:
            exposure.getPsf().setCacheSize(psfCache)

        allMeasCat = None

        for tract, refCat in tractRefCat.items():
            if refCat is None:
                self.log.warn('Failed to get references for %s. skipping' % dataRef.dataId)
                continue
            refWcs = tract.getWcs()

            measCat = self.measurement.generateMeasCat(exposure, refCat, refWcs,
                                                       idFactory=self.makeIdFactory(dataRef))

            self.log.info("Performing forced measurement on %s" % (dataRef.dataId,))
            self.attachFootprints(measCat, refCat, exposure, refWcs, dataRef)

            self.measurement.run(measCat, exposure, refCat, refWcs, exposureId=self.getExposureId(dataRef))

            if self.config.doApCorr:
                self.applyApCorr.run(
                    catalog=measCat,
                    apCorrMap=exposure.getInfo().getApCorrMap()
                )
            self.catalogCalculation.run(measCat)

            if allMeasCat is None:
                allMeasCat = measCat
            else:
                allMeasCat.extend(measCat)

        if allMeasCat is not None:
            self.writeOutput(dataRef, allMeasCat)

    def getReferences(self, butler, exposure):
        """!
        Get reference catalogs from all patch,tract combinations that overlaps this exposure

        This method will filter out any rferences that are not in the inner tract,patch region
        so that there will be no duplicates across boundaries.
        @param[in]  butler     An lsst.daf.persistence.Butler.  This is used to get
                               the references
        @param[in]  exposure   A deepDiff_exposure on which to run the measurements

        @return    Dictionary of tract, reference catalogs

        """
        skyMap = butler.get(f"{self.config.coaddName}Coadd_skyMap")

        wcs = exposure.getWcs()
        imagePixelCorners = exposure.getBBox().getCorners()
        imageSkyCorners = [wcs.pixelToSky(afwGeom.Point2D(a)) for a in imagePixelCorners]

        # get initial list of tracts and patches
        tractList = skyMap.findTractPatchList(imageSkyCorners)

        tractRefCat = dict()
        for tract, patchList in tractList:
            usePatchList = []
            for patch in patchList:
                patchPoly = patch.getInnerSkyPolygon(tract.getWcs())
                imagePoly = sphgeom.ConvexPolygon.convexHull([coord.getVector() for coord in imageSkyCorners])
                if patchPoly.intersects(imagePoly):
                    usePatchList.append(patch)

            refCat = self.fetchInPatches(butler, exposure, tract, usePatchList)
            tractRefCat[tract] = refCat

        return tractRefCat

    def fetchInPatches(self, butler, exposure, tract, patchList):
        """!
        Get the reference catalogs from a given tract,patchlist

        @param[in]  butler     A Butler used to get the reference catalogs
        @param[in]  exposure   A deepDiff_exposure on which to run the measurements
        @param[in]  tract      The tract
        @param[in]  patchList  A list of patches that need to be checked


        @return    Combined SourceCatalog from all the patches
        """
        dataset = f"{self.config.coaddName}Diff_diaObject"
        catalog = None

        for patch in patchList:
            dataId = {'tract': tract.getId(), 'patch': "%d,%d" % patch.getIndex()}
            self.log.info("Getting references in %s" % (dataId,))
            if not butler.datasetExists(dataset, dataId):
                if self.config.skipMissing:
                    self.log.info("Could not find %s for dataset %s" % (dataId, dataset))
                    continue
                raise pipeBase.TaskError("Reference %s doesn't exist" % (dataId,))

            new_catalog = butler.get(dataset, dataId, immediate=True)
            patchBox = geom.Box2D(patch.getInnerBBox())
            tractBox = tract.getInnerSkyPolygon()
            tractWcs = tract.getWcs()
            expBox = geom.Box2D(exposure.getBBox())
            expWcs = exposure.getWcs()

            # only use objects that overlap inner patch bounding box and overlap exposure
            validPatch = np.array([
                patchBox.contains(tractWcs.skyToPixel(s.getCoord())) for s in new_catalog])

            # There doesn't seem to be a inner bounding box so I have to use the sphgeom stuff
            validTract = np.array(
                [tractBox.contains(sphgeom.UnitVector3d(sphgeom.LonLat.fromRadians(s.getRa().asRadians(),
                                                                                   s.getDec().asRadians())))
                    for s in new_catalog])
            validExposure = np.array([
                expBox.contains(expWcs.skyToPixel(s.getCoord())) for s in new_catalog])

            if validPatch.size == 0 or validExposure.size == 0 or validTract.size == 0:
                self.log.debug("No valid sources %s for dataset %s" % (dataId, dataset))
                continue

            if catalog is None:
                catalog = new_catalog[validPatch & validExposure & validTract]
            else:
                catalog.extend(new_catalog[validPatch & validExposure & validTract])

        return catalog


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

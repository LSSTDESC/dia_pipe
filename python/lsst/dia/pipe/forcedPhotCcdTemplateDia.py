import numpy as np

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom as geom
import lsst.sphgeom as sphgeom

from lsst.meas.base.forcedPhotCcd import ForcedPhotCcdTask, ForcedPhotCcdConfig

from .forcedPhotDia import DiaReferencesTask


__all__ = ("ForcedPhotCcdTemplateDiaConfig", "ForcedPhotCcdTemplateDiaTask")


class ForcedPhotCcdTemplateDiaConfig(ForcedPhotCcdConfig):
    skipMissing = pexConfig.Field(dtype=bool, default=True,
                                  doc="Skip getting references if they do not exist?")
    coaddName = pexConfig.Field(dtype=str, default='deep',
                                doc="Name of coadd")

    def setDefaults(self):
        ForcedPhotCcdTask.ConfigClass.setDefaults(self)
        self.references.retarget(DiaReferencesTask)
        self.references.datasetSuffix = 'deepCoadd_meas'

        self.measurement.copyColumns = {"id": "template_id", "coord_ra": "coord_ra",
                                        "coord_dec": "coord_dec",
                                        "base_PsfFlux_instFlux":"template_base_PsfFlux_instFlux",
                                        "base_PsfFlux_instFluxErr": "template_base_PsfFlux_instFluxErr",
                                        "detect_isPrimary": "detect_isPrimary"}

        self.measurement.plugins.names = ['base_SdssShape', 'base_DiaTransformedCentroid',
                                          'base_PsfFlux', 'base_LocalBackground',
                                          'base_PixelFlags', 'base_CircularApertureFlux']
        self.measurement.slots.centroid = 'base_DiaTransformedCentroid'
        self.measurement.slots.shape = 'base_SdssShape'
        self.measurement.slots.apFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.psfFlux = None
        self.measurement.slots.calibFlux = None

        # These radii were chosen because they are among the default measured in the pipeline.  If the default
        # changes then these will not be able to be copied.
        radii = [3., 6., 9., 12.]
        for radius in radii:
            base = int(radius)
            decimal = int((radius - int(radius))*10)
            input_name = f"base_CircularApertureFlux_{base}_{decimal}_instFlux"
            output_name = f"template_base_CircularApertureFlux_{base}_{decimal}_instFlux"
            self.measurement.copyColumns[input_name] = output_name

            input_name = f"base_CircularApertureFlux_{base}_{decimal}_instFluxErr"
            output_name = f"template_base_CircularApertureFlux_{base}_{decimal}_instFluxErr"
            self.measurement.copyColumns[input_name] = output_name

        self.measurement.plugins["base_CircularApertureFlux"].radii = radii
        # Use a large aperture to be independent of seeing in calibration
        self.measurement.plugins["base_CircularApertureFlux"].maxSincRadius = 12.0


class ForcedPhotCcdTemplateDiaTask(ForcedPhotCcdTask):
    """!
    A command-line driver for performing forced measurement on difference images from template catalogs.

    It will look for all tract/patches that overlap the bounding box.  If the parent of an object
    does not fall within the bounding box, then the object will be removed.
    """

    ConfigClass = ForcedPhotCcdTemplateDiaConfig
    RunnerClass = pipeBase.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCcdTemplateDia"
    dataPrefix = "deepDiff_"


    def __init__(self, butler=None, refSchema=None, **kwds):
        """Initialize the task.
        ForcedPhotImageTask takes two keyword arguments beyond the usual CmdLineTask arguments:
         - refSchema: the Schema of the reference catalog, passed to the constructor of the references
           subtask
         - butler: a butler that will be passed to the references subtask to allow it to load its Schema
           from disk
        At least one of these arguments must be present; if both are, schema takes precedence.
        """
        super(ForcedPhotCcdTask, self).__init__(butler, refSchema, **kwds)

    def writeOutput(self, dataRef, sources):
        """!Write source table
        @param dataRef  Data reference from butler
        @param sources  SourceCatalog to save
        """
        dataRef.put(sources, self.dataPrefix + "forced_diaSrcTemplate",
                    flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

    def _getConfigName(self):
        """!Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return "forcedPhotCcdTemplateDia_config"

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
        tractRefCat = self.getReferences(dataRef.getButler(), exposure, dataRef.dataId['filter'])

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

    def getReferences(self, butler, exposure, band):
        """!
        Get reference catalogs from all patch,tract combinations that overlaps this exposure

        This method will filter out any rferences that are not in the outer patch region
        so there will be duplicates across boundaries.  These can be filtered out
        using the detect_isPrimary flag.

        @param[in]  butler     An lsst.daf.persistence.Butler.  This is used to get
                               the references
        @param[in]  exposure   A deepDiff_exposure on which to run the measurements

        @return    Dictionary of tract, reference catalogs

        """
        skyMap = butler.get(f"{self.config.coaddName}Coadd_skyMap")

        wcs = exposure.getWcs()
        imagePixelCorners = exposure.getBBox().getCorners()
        imageSkyCorners = [wcs.pixelToSky(geom.Point2D(a)) for a in imagePixelCorners]

        # get initial list of tracts and patches
        tractList = skyMap.findTractPatchList(imageSkyCorners)

        tractRefCat = dict()
        for tract, patchList in tractList:
            usePatchList = []
            for patch in patchList:
                patchPoly = patch.getOuterSkyPolygon(tract.getWcs())
                imagePoly = sphgeom.ConvexPolygon.convexHull([coord.getVector() for coord in imageSkyCorners])
                if patchPoly.intersects(imagePoly):
                    usePatchList.append(patch)

            refCat = self.fetchInPatches(butler, exposure, tract, usePatchList, band)
            if refCat:
                refCat.sort(afwTable.SourceTable.getParentKey())
            tractRefCat[tract] = refCat

        return tractRefCat

    def fetchInPatches(self, butler, exposure, tract, patchList, band):
        """!
        Get the reference catalogs from a given tract,patchlist
        This will remove objects where the child is inside the catalog boundary, but
        the parent is outside the boundary.

        @param[in]  butler     A Butler used to get the reference catalogs
        @param[in]  exposure   A deepDiff_exposure on which to run the measurements
        @param[in]  tract      The tract
        @param[in]  patchList  A list of patches that need to be checked
.

        @return    Combined SourceCatalog from all the patches
        """
        dataset = f"{self.config.coaddName}Coadd_meas"
        catalog = None

        for patch in patchList:
            dataId = {'tract': tract.getId(), 'patch': "%d,%d" % patch.getIndex()}
                    
            
            dataId['filter'] = band
            self.log.info("Getting references in %s" % (dataId,))
            if not butler.datasetExists(dataset, dataId):
                if self.config.skipMissing:
                    self.log.info("Could not find %s for dataset %s" % (dataId, dataset))
                    continue
                raise pipeBase.TaskError("Reference %s doesn't exist" % (dataId,))

            new_catalog = butler.get(dataset, dataId, immediate=True)
            patchBox = geom.Box2D(patch.getOuterBBox())
            tractWcs = tract.getWcs()
            expBox = geom.Box2D(exposure.getBBox())
            expWcs = exposure.getWcs()

            # only use objects that overlap patch bounding box and overlap exposure
            validPatch = np.array([
                patchBox.contains(tractWcs.skyToPixel(s.getCoord())) for s in new_catalog])

            validExposure = np.array([
                expBox.contains(expWcs.skyToPixel(s.getCoord())) for s in new_catalog])

            if validPatch.size == 0 or validExposure.size == 0:
                self.log.debug("No valid sources %s for dataset %s" % (dataId, dataset))
                continue

            if catalog is None:
                catalog = new_catalog[validPatch & validExposure]
            else:
                catalog.extend(new_catalog[validPatch & validExposure])

        if catalog is None:
            return None
        # if the parent is not in the catalog remove
        refCatIdDict = {ref.getId(): ref.getParent() for ref in catalog}
        refCatIdDict[0] = 0
        parentGood = np.array([refCatIdDict[ref.getId()] in refCatIdDict  for ref in catalog])
        if np.sum(parentGood==False) > 1:
            self.log.info("Removing %d/%d objects without parents" % (np.sum(parentGood==False),
                len(parentGood)))
            catalog = catalog.copy(deep=True)[parentGood]

        return catalog

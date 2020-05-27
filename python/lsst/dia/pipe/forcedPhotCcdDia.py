import numpy as np

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom as geom
import lsst.sphgeom as sphgeom

from lsst.meas.base.forcedPhotCcd import ForcedPhotCcdTask, ForcedPhotCcdConfig

from .forcedPhotDia import DiaReferencesTask

__all__ = ("ForcedPhotCcdDiaConfig", "ForcedPhotCcdDiaTask")


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
                              references subtask to obtain the reference WCS, the image, 
                              and load the reference catalog.    Sources are
                              generated with generateMeasCat() in the measurement subtask.  These
                              are passed to measurement's run method which fills the source
                              catalog with the forced measurement results.  The sources are then
                              passed to the writeOutputs() method (implemented by derived classes)
                              which writes the outputs.  See derived class documentation for which
                              datasets and data ID keys are used.
        @param[in]  psfCache  Size of PSF cache, or None. The size of the PSF cache can have
                              a significant effect upon the runtime for complicated PSF models.
        """
        exposure = dataRef.get('deepDiff_differenceExp')
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
        imageSkyCorners = [wcs.pixelToSky(geom.Point2D(a)) for a in imagePixelCorners]

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


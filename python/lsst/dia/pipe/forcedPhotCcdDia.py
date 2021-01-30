import numpy as np
import os

import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
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


class ForcedPhotCatalogDiaConfig(ForcedPhotCcdConfig):
    footprintSize = pexConfig.Field(dtype=float, default=10,
                                doc=" multiply this by the size of the PSF to determine the footprint")
    catalogName = pexConfig.Field(dtype=str, default='cat.fits',
                                doc="Name of input catalog of forced measurements")
    outputDir = pexConfig.Field(dtype=str, default='.',
                                doc="Name of output directory")
    outputName = pexConfig.Field(dtype=str, default='forced_{visit}_{detector}.fits',
                                doc="Name of output catalog of forced measurements")

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


class ForcedPhotCatalogDiaTask(ForcedPhotCcdTask):
    """!
    A command-line driver for performing forced measurement on CCD images from a catalog of object locations.
    The catalog should contain three columns 'coord_ra' and 'coord_dec' in radians and an id.

    This task is a subclass of ForcedPhotCcdTask, although it does not use most of the
    functionality defined there.  In this task we delegate the process of looking up
    the reference catalog inside the run function.  We allow for the references to come
    from multiple tracts and perform forced phometry on each tract and then combine them.
    The refernces are all taken from the inner patch boundary
    """

    ConfigClass = ForcedPhotCatalogDiaConfig
    RunnerClass = pipeBase.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCatalogDia"
    dataPrefix = ""

    def writeOutput(self, dataRef, sources):
        """!Write source table
        @param dataRef  Data reference from butler
        @param sources  SourceCatalog to save
        """

        name = self.config.outputName.format(**dataRef.dataId)
        sources.writeFits(os.path.join(self.config.outputDir, name))

    def _getConfigName(self):
        """!Return the name of the config dataset.
        """
        return None

    def _getMetadataName(self):
        """!Return the name of the metadata dataset.
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
        psfSize = exposure.getPsf().computeShape().getDeterminantRadius()
        refWcs = exposure.getWcs()

        refCatBase = afwTable.BaseCatalog.readFits(self.config.catalogName)
        refSchema = afwTable.SourceTable.makeMinimalSchema()
        refCat = afwTable.SourceCatalog(refSchema)
        for ref in refCatBase:
            nref = refCat.addNew()
            xy = refWcs.skyToPixel(geom.SpherePoint(ref['coord_ra']*geom.radians, ref['coord_dec']*geom.radians))
            bbox = geom.Box2I(geom.Point2I(int(xy.getX()), int(xy.getY())), geom.Extent2I(1, 1))
            bbox.grow(int(self.config.footprintSize*psfSize))
            footprint = afwDet.Footprint(afwGeom.SpanSet(bbox), bbox)
            nref.setFootprint(footprint)
            nref.setParent(0)
            nref.set('coord_ra', ref['coord_ra']*geom.radians)
            nref.set('coord_dec', ref['coord_dec']*geom.radians)
            nref.set('id', ref['id'])
            

        if psfCache is not None:
            exposure.getPsf().setCacheSize(psfCache)

        measCat = self.generateMeasCat(exposure, refCat, refWcs,
                                       idFactory=self.makeIdFactory(dataRef))

        self.log.info("Performing forced measurement on %s" % (dataRef.dataId,))
        self.attachFootprints(measCat, refCat, exposure, refWcs, dataRef)

        self.run(measCat, exposure, refCat, refWcs, exposureId=self.getExposureId(dataRef))

        if self.config.doApCorr:
            self.applyApCorr.run(
                catalog=measCat,
                apCorrMap=exposure.getInfo().getApCorrMap()
            )
        self.catalogCalculation.run(measCat)

        if measCat is not None:
            self.writeOutput(dataRef, measCat)

    def generateMeasCat(self, exposure, refCat, refWcs, idFactory=None):
        """!Initialize an output SourceCatalog using information from the reference catalog.

        Since are using a potentially different type of input catalog we overide the
        default catalog which uses a schema mapper object to copy information.

        This generates a new blank SourceRecord for each record in refCat.  Note that this
        method does not attach any Footprints.  Doing so is up to the caller (who may
        call attachedTransformedFootprints or define their own method - see run() for more
        information).
        @param[in] exposure    Exposure to be measured
        @param[in] refCat      Sequence (not necessarily a SourceCatalog) of reference SourceRecords.
        @param[in] refWcs      Wcs that defines the X,Y coordinate system of refCat
        @param[in] idFactory   factory for creating IDs for sources
        @return    Source catalog ready for measurement
        """
        if idFactory is None:
            idFactory = afwTable.IdFactory.makeSimple()
        table = afwTable.SourceTable.make(self.measurement.schema, idFactory)
        measCat = afwTable.SourceCatalog(table)
        table = measCat.table
        table.setMetadata(self.measurement.algMetadata)
        table.preallocate(len(refCat))

        for ref in refCat:
            newSource = measCat.addNew()
            for key,val in  self.config.measurement.copyColumns.items():
                newSource.set(val, ref[key])
            
        return measCat
        
    def attachFootprints(self, sources, refCat, exposure, refWcs, dataRef):
        """Copy footprints from reference catalog to output catalog
        """

        for src,ref in zip(sources, refCat):
            src.setFootprint(ref.getFootprint())
    
    def run(self, measCat, exposure, refCat, refWcs, exposureId=None, beginOrder=None, endOrder=None):
        """!
        Perform forced measurement.
        Simplified measurement that does not try to insert blended objects

        @param[in]  exposure     lsst.afw.image.ExposureF to be measured; must have at least a Wcs attached.
        @param[in]  measCat      Source catalog for measurement results; must be initialized with empty
                                 records already corresponding to those in refCat (via e.g. generateMeasCat).
        @param[in]  refCat       A sequence of SourceRecord objects that provide reference information
                                 for the measurement.  These will be passed to each Plugin in addition
                                 to the output SourceRecord.
        @param[in]  refWcs       Wcs that defines the X,Y coordinate system of refCat
        @param[in]  exposureId   optional unique exposureId used to calculate random number
                                 generator seed in the NoiseReplacer.
        @param[in]  beginOrder   beginning execution order (inclusive): measurements with
                                 executionOrder < beginOrder are not executed. None for no limit.
        @param[in]  endOrder     ending execution order (exclusive): measurements with
                                 executionOrder >= endOrder are not executed. None for no limit.
        Fills the initial empty SourceCatalog with forced measurement results.  Two steps must occur
        before run() can be called:
         - generateMeasCat() must be called to create the output measCat argument.
         - Footprints appropriate for the forced sources must be attached to the measCat records.  The
           attachTransformedFootprints() method can be used to do this, but this degrades HeavyFootprints
           to regular Footprints, leading to non-deblended measurement, so most callers should provide
           Footprints some other way.  Typically, calling code will have access to information that will
           allow them to provide HeavyFootprints - for instance, ForcedPhotCoaddTask uses the HeavyFootprints
           from deblending run in the same band just before non-forced is run measurement in that band.
        """

        for ii,(refRecord, measRecord) in enumerate(zip(refCat, measCat)):
            self.measurement.callMeasure(measRecord, exposure, refRecord, refWcs,
                                        beginOrder=beginOrder, endOrder=endOrder)


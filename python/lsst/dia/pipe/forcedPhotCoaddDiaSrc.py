import numpy as np

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom as geom
import lsst.sphgeom as sphgeom

from lsst.meas.base.forcedPhotCcd import ForcedPhotCcdTask, ForcedPhotCcdConfig
from .forcedPhotDia import DiaSrcReferencesTask

__all__ = ("ForcedPhotCoaddDiaSrcConfig", "ForcedPhotCoaddDiaSrcTask")

class ForcedPhotCoaddDiaSrcConfig(ForcedPhotCcdConfig):
    coaddName = pexConfig.Field(dtype=str, default='deep',
                                doc="Name of coadd")

    def setDefaults(self):
        ForcedPhotCcdTask.ConfigClass.setDefaults(self)
        self.references.retarget(DiaSrcReferencesTask)
        self.measurement.copyColumns = {"id": "id", "coord_ra": "coord_ra",
                                        "coord_dec": "coord_dec",
                                        "base_PsfFlux_instFlux":"diaSrc_base_PsfFlux_instFlux",
                                        "base_PsfFlux_instFluxErr": "diaSrc_base_PsfFlux_instFluxErr",
                                        }
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
            output_name = f"diaSrc_base_CircularApertureFlux_{base}_{decimal}_instFlux"
            self.measurement.copyColumns[input_name] = output_name

            input_name = f"base_CircularApertureFlux_{base}_{decimal}_instFluxErr"
            output_name = f"diaSrc_base_CircularApertureFlux_{base}_{decimal}_instFluxErr"
            self.measurement.copyColumns[input_name] = output_name

        self.measurement.plugins["base_CircularApertureFlux"].radii = radii
        # Use a large aperture to be independent of seeing in calibration
        self.measurement.plugins["base_CircularApertureFlux"].maxSincRadius = 12.0


class ForcedPhotCoaddDiaSrcTask(ForcedPhotCcdTask):
    """!
    A command-line driver for performing forced measurement on Coadd images from DIASrc catalogs.


    """

    ConfigClass = ForcedPhotCoaddDiaSrcConfig
    RunnerClass = pipeBase.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCoaddDiaSrc"
    dataPrefix = "deepCoadd_"

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
        self.primaryKey = self.measurement.schema.addField("detect_isPrimary", type="Flag", doc="set to True if inside inner patch and tract region")

    def writeOutput(self, dataRef, sources):
        """!Write source table
        @param dataRef  Data reference from butler
        @param sources  SourceCatalog to save
        """
        dataRef.put(sources, "deepDiff_forced_template_diaSrc",
                    flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

    def _getConfigName(self):
        """!Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return "forcedPhotCoaddDiaSrc_config"

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

    def runDataRef(self, dataRef):
        """!Measure a single exposure using forced detection for a reference catalog.
        @param[in]  dataRef   An lsst.daf.persistence.ButlerDataRef
        @param[in]  psfCache  Size of PSF cache, or None. The size of the PSF cache can have
                              a significant effect upon the runtime for complicated PSF models.
        """
        exposure = dataRef.get('deepDiff_differenceExp')
        catalog = dataRef.get('deepDiff_diaSrc')
        expWcs = exposure.getWcs()
        butler = dataRef.butlerSubset.butler

        # I need to get the template images/catalogs for all overlapping tracts
        skyMap = butler.get(datasetType=self.config.coaddName + "Coadd_skyMap")
        skyCorners = [expWcs.pixelToSky(geom.Point2D(pixPos)) for pixPos in exposure.getBBox().getCorners()]
        imagePoly = sphgeom.ConvexPolygon.convexHull([coord.getVector() for coord in skyCorners])
        tractPatchList = skyMap.findTractPatchList(skyCorners)

        allMeasCat = None
        for tract, patchList in tractPatchList:
            for patch in patchList:
                self.log.info('Processing patch %s from tract %s' % (patch.getIndex(),tract))
                patchPoly = patch.getOuterSkyPolygon(tract.getWcs())
                if patchPoly.intersects(imagePoly) is False:
                    self.log.info('No intersection with the boundary patch.')
                    continue

                validObject = np.array(
                [patchPoly.contains(sphgeom.UnitVector3d(sphgeom.LonLat.fromRadians(s.getRa().asRadians(),
                                                                                    s.getDec().asRadians())))
                 for s in catalog])

                refCat = catalog[validObject]
                expCorners = [tract.getWcs().skyToPixel(pos) for pos in skyCorners]
                expBBox = geom.Box2I()
                for corner in expCorners:
                    expBBox.include(geom.Point2I(corner))

                overlapBox = geom.Box2I(patch.getOuterBBox())
                overlapBox.clip(expBBox)
                
                patchArgDict = dict(
                        datasetType=self.dataPrefix+ "calexp_sub",
                        bbox=overlapBox,
                        tract=tract.getId(),
                        patch="%s,%s" % (patch.getIndex()[0], patch.getIndex()[1]),
                        filter=exposure.getFilter().getName()
                    )

                coaddPatch = butler.get(**patchArgDict)

                # I need to filter out objects whose parents are not in the bounding box
                refCatIdDict = {ref.getId(): ref.getParent() for ref in refCat}
                # Add 0 for objects without a parent
                refCatIdDict[0] = 0    
                parentGood = np.array([refCatIdDict[ref.getId()] in refCatIdDict  for ref in refCat])

                if np.sum(parentGood==False) > 1:
                    self.log.info("Removing %d/%d objects without parents" % (np.sum(parentGood==False),
                                  len(parentGood)))
                    refCat = refCat.copy(deep=True)[parentGood]

                if len(refCat) == 0:
                    self.log.info('No references available.')
                    continue

                measCat = self.measurement.generateMeasCat(coaddPatch, refCat, expWcs,
                                                           idFactory=self.makeIdFactory(dataRef))

                self.log.info("Performing forced measurement on %s" % (patchArgDict))
                self.attachFootprints(measCat, refCat, coaddPatch, expWcs, dataRef)

                exposureId = self.getExposureId(dataRef)
                self.measurement.run(measCat, coaddPatch, refCat, expWcs, exposureId=exposureId)

                # Label primary objects
                innerBox = geom.Box2D(patch.getInnerBBox())
                insideBox = np.array([innerBox.contains(s.getCentroid()) for s in measCat])
                primaryTract = np.array([skyMap.findTract(s.getCoord())==tract for s in measCat])
                primary = (insideBox) & (primaryTract)
                # I can't set the whole array, so I do it one item at a time
                for s,p in zip(measCat, primary):
                    s.set(self.primaryKey,p)

                if self.config.doApCorr:
                    self.applyApCorr.run(
                        catalog=measCat,
                        apCorrMap=coaddPatch.getInfo().getApCorrMap()
                    )
                self.catalogCalculation.run(measCat)

                if allMeasCat is None:
                    allMeasCat = measCat
                else:
                    allMeasCat.extend(measCat)

        if allMeasCat is not None:
            self.writeOutput(dataRef, allMeasCat)

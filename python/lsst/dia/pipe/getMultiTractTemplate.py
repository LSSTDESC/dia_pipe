import numpy as np

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom

from lsst.meas.algorithms import CoaddPsf, CoaddPsfConfig
from lsst.meas.algorithms import WarpedPsf


__all__ = ["GetCoaddAsMultiTractTemplateTask", "GetCoaddAsMultiTractTemplateConfig"]


class GetCoaddAsMultiTractTemplateConfig(pexConfig.Config):
    templateBorderSize = pexConfig.Field(
        dtype=int,
        default=20,
        doc="Number of pixels to grow the requested template image to account for warping"
    )
    coaddName = pexConfig.Field(
        doc="coadd name: typically one of 'deep', 'goodSeeing', or 'dcr'",
        dtype=str,
        default="deep",
    )
    warpType = pexConfig.Field(
        doc="Warp type of the coadd template: one of 'direct' or 'psfMatched'",
        dtype=str,
        default="direct",
    )
    coaddPsf = pexConfig.ConfigField(
        doc="Configuration for CoaddPsf",
        dtype=CoaddPsfConfig,
    )
    warp = pexConfig.ConfigField(
        dtype=afwMath.Warper.ConfigClass,
        doc="warper configuration",
    )
    statistic = pexConfig.Field(
        dtype=str,
        doc="How to combine tracts that overlap",
        default="MEAN",
    )


class GetCoaddAsMultiTractTemplateTask(pipeBase.Task):
    """Subtask to retrieve coadd from possibly different tracts and 
    use as an image difference template.  It uses the tract closest to the
    central point of the ccd as the reference tract.  All other tracts will
    be warped onto the reference task.

    The PSF of the resulting template will be a CoaddPSF of individual CoaddPSFs.
    """

    ConfigClass = GetCoaddAsMultiTractTemplateConfig
    _DefaultName = "GetCoaddAsMultiTractTemplateTask"

    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.warper = afwMath.Warper.fromConfig(self.config.warp)


    def run(self, exposure, sensorRef, templateIdList=None):
        """Retrieve and mosaic a template coadd that overlaps the exposure where
        the template spans multiple tracts.

        The resulting template image will be an average of all the input templates from
        the separate tracts.

        The PSF on the template is created by combining the CoaddPsf on each template image
        into a meta-CoaddPsf.

        Parameters
        ----------
        exposure: `lsst.afw.image.Exposure`
            an exposure for which to generate an overlapping template
        sensorRef : TYPE
            a Butler data reference that can be used to obtain coadd data
        templateIdList : TYPE, optional
            list of data ids (unused)
        Returns
        -------
        result : `struct`
            return a pipeBase.Struct:
            - ``exposure`` : a template coadd exposure assembled out of patches
            - ``sources`` :  None for this subtask
        """

        # Table for CoaddPSF
        tractsSchema = afwTable.ExposureTable.makeMinimalSchema()
        tractKey = tractsSchema.addField('tract', type=np.int32, doc='Which tract')
        patchKey = tractsSchema.addField('patch', type=np.int32, doc='Which patch')
        weightKey = tractsSchema.addField('weight', type=float, doc='Weight for each tract, should be 1')
        tractsCatalog = afwTable.ExposureCatalog(tractsSchema)


        skyMap = sensorRef.get(datasetType=self.config.coaddName + "Coadd_skyMap")
        expWcs = exposure.getWcs()
        expBoxD = geom.Box2D(exposure.getBBox())
        expBoxD.grow(self.config.templateBorderSize)
        ctrSkyPos = expWcs.pixelToSky(expBoxD.getCenter())

        centralTractInfo = skyMap.findTract(ctrSkyPos)
        if not centralTractInfo:
            raise RuntimeError("No suitable tract found for central point")

        self.log.info("Central skyMap tract %s" % (centralTractInfo.getId(),))

        skyCorners = [expWcs.pixelToSky(pixPos) for pixPos in expBoxD.getCorners()]
        tractPatchList = skyMap.findTractPatchList(skyCorners)
        if not tractPatchList:
            raise RuntimeError("No suitable tract found")

        self.log.info("All overlapping skyMap tracts %s" % ([a[0].getId() for a in tractPatchList]))

        # Move central tract to front of the list and use as the reference
        tracts = [tract[0].getId() for tract in tractPatchList]
        centralIndex = tracts.index(centralTractInfo.getId())
        tracts.insert(0, tracts.pop(centralIndex))
        tractPatchList.insert(0, tractPatchList.pop(centralIndex))

        coaddPsf = None
        coaddFilter = None
        nPatchesFound = 0

        maskedImageList=[]
        weightList=[]

        for itract,tract in enumerate(tracts):
            tractInfo = tractPatchList[itract][0]

            coaddWcs = tractInfo.getWcs()
            coaddBBox = geom.Box2D()
            for skyPos in skyCorners:
                coaddBBox.include(coaddWcs.skyToPixel(skyPos))
            coaddBBox = geom.Box2I(coaddBBox)
                
            if itract == 0:
                # Define final wcs and bounding box from the reference tract
                finalWcs = coaddWcs
                finalBBox = coaddBBox    

            patchList = tractPatchList[itract][1]
            for patchInfo in patchList:
                self.log.info('Adding patch %s from tract %s' % (patchInfo.getIndex(),tract))

                # Local patch information
                patchSubBBox = geom.Box2I(patchInfo.getInnerBBox())
                patchSubBBox.clip(coaddBBox)
                patchInt = int(f"{patchInfo.getIndex()[0]}{patchInfo.getIndex()[1]}")
                innerBBox = geom.Box2I(tractInfo._minimumBoundingBox(finalWcs))

                if itract == 0:          
                    # clip to image and tract boundaries  
                    patchSubBBox.clip(finalBBox)
                    patchSubBBox.clip(innerBBox)
                    if patchSubBBox.getArea() == 0:
                        self.log.debug("No ovlerap for patch %s" % patchInfo)
                        continue

                    patchArgDict = dict(
                        datasetType="deepCoadd_sub",
                        bbox=patchSubBBox,
                        tract=tractInfo.getId(),
                        patch="%s,%s" % (patchInfo.getIndex()[0], patchInfo.getIndex()[1]),
                        filter=exposure.getFilter().getName()
                    )
                    coaddPatch = sensorRef.get(**patchArgDict)
                    if coaddFilter is None:
                        coaddFilter = coaddPatch.getFilter()

                    # create full image from final bounding box
                    exp = afwImage.ExposureF(finalBBox, finalWcs)
                    exp.maskedImage.set(np.nan, afwImage.Mask.getPlaneBitMask("NO_DATA"), np.nan)
                    exp.maskedImage.assign(coaddPatch.maskedImage, patchSubBBox)

                    maskedImageList.append(exp.maskedImage)
                    weightList.append(1)

                    record = tractsCatalog.addNew()
                    record.setPsf(coaddPatch.getPsf())
                    record.setWcs(coaddPatch.getWcs())
                    record.setPhotoCalib(coaddPatch.getPhotoCalib())
                    record.setBBox(patchSubBBox)
                    record.set(tractKey, tract)
                    record.set(patchKey, patchInt)
                    record.set(weightKey, 1.)
                    nPatchesFound += 1
                else:
                    # compute the exposure bounding box in a tract that is not the reference tract
                    localBox = geom.Box2I()
                    for skyPos in skyCorners:
                        localBox.include(geom.Point2I(tractInfo.getWcs().skyToPixel(skyPos)))
                    
                    # clip to patch bounding box
                    localBox.clip(patchSubBBox)

                    # grow border to deal with warping at edges
                    localBox.grow(self.config.templateBorderSize)

                    # clip to tract inner bounding box
                    localInnerBBox = geom.Box2I(tractInfo._minimumBoundingBox(tractInfo.getWcs()))
                    localBox.clip(localInnerBBox)

                    patchArgDict = dict(
                        datasetType="deepCoadd_sub",
                        bbox=localBox,
                        tract=tractInfo.getId(),
                        patch="%s,%s" % (patchInfo.getIndex()[0], patchInfo.getIndex()[1]),
                        filter=exposure.getFilter().getName()
                    )
                    coaddPatch = sensorRef.get(**patchArgDict)

                    # warp to reference tract wcs
                    xyTransform = afwGeom.makeWcsPairTransform(coaddPatch.getWcs(), finalWcs)
                    psfWarped = WarpedPsf(coaddPatch.getPsf(), xyTransform)
                    warped = self.warper.warpExposure(finalWcs, coaddPatch, maxBBox=finalBBox)

                    # check if warpped image is viable
                    if warped.getBBox().getArea() == 0:
                        self.log.info("No ovlerap for warped patch %s. Skipping" % patchInfo)
                        continue

                    warped.setPsf(psfWarped)

                    exp = afwImage.ExposureF(finalBBox, finalWcs)
                    exp.maskedImage.set(np.nan, afwImage.Mask.getPlaneBitMask("NO_DATA"), np.nan)
                    exp.maskedImage.assign(warped.maskedImage, warped.getBBox())

                    maskedImageList.append(exp.maskedImage)
                    weightList.append(1)
                    record = tractsCatalog.addNew()
                    record.setPsf(psfWarped)
                    record.setWcs(finalWcs)
                    record.setPhotoCalib(coaddPatch.getPhotoCalib())
                    record.setBBox(warped.getBBox())
                    record.set(tractKey, tract)
                    record.set(patchKey, patchInt)
                    record.set(weightKey, 1.)
                    nPatchesFound += 1

         
        if nPatchesFound == 0:
            raise RuntimeError("No patches found!")

        coaddPsf = CoaddPsf(tractsCatalog, coaddWcs, self.config.coaddPsf.makeControl())

        # Combine images from individual patches together

        # Do not mask any values
        statsFlags = afwMath.stringToStatisticsProperty(self.config.statistic)
        maskMap = []
        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNanSafe(True)
        statsCtrl.setWeighted(True)
        statsCtrl.setCalcErrorFromInputVariance(True)
        
        coaddExposure = afwImage.ExposureF(finalBBox, finalWcs)
        coaddExposure.maskedImage.set(np.nan, afwImage.Mask.getPlaneBitMask("NO_DATA"), np.nan)
        xy0 = coaddExposure.getXY0()
        coaddExposure.maskedImage = afwMath.statisticsStack(maskedImageList, 
            statsFlags, statsCtrl, weightList, 0, maskMap)
        coaddExposure.maskedImage.setXY0(xy0)   
        coaddExposure.writeFits('test.fits')
        coaddPsf = CoaddPsf(tractsCatalog, finalWcs, self.config.coaddPsf.makeControl())
        if coaddPsf is None:
            raise RuntimeError("No coadd Psf found!")

        coaddExposure.setPsf(coaddPsf)
        coaddExposure.setFilter(coaddFilter)
        return pipeBase.Struct(exposure=coaddExposure, sources=None)

    def getCoaddDatasetName(self):
        """Return coadd name for given task config
        Returns
        -------
        CoaddDatasetName : `string`
        TODO: This nearly duplicates a method in CoaddBaseTask (DM-11985)
        """
        warpType = self.config.warpType
        suffix = "" if warpType == "direct" else warpType[0].upper() + warpType[1:]
        return self.config.coaddName + "Coadd" + suffix



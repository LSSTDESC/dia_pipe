import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.daf.base import DateTime
from lsst.pipe.tasks.selectImages import WcsSelectImagesTask

__all__ = ["TimeBestSeeingWcsSelectImagesTask"]


class TimeBestSeeingWcsSelectImageConfig(WcsSelectImagesTask.ConfigClass):
    """Base configuration for BestSeeingSelectImagesTask.
    """
    nImagesMax = pexConfig.Field(
        dtype=int,
        doc="Maximum number of images to select, ignored if None",
        default=None)
    maxPsfFwhm = pexConfig.Field(
        dtype=float,
        doc="Maximum PSF FWHM (in pixels) to select",
        default=None,
        optional=True)
    minPsfFwhm = pexConfig.Field(
        dtype=float,
        doc="Minimum PSF FWHM (in pixels) to select",
        default=None,
        optional=True)
    minMJD = pexConfig.Field(
        dtype=float,
        doc="Minimum MJD to select",
        default=None,
        optional=True)
    maxMJD = pexConfig.Field(
        dtype=float,
        doc="Maximum MJD to select",
        default=None,
        optional=True)


class TimeBestSeeingWcsSelectImagesTask(WcsSelectImagesTask):
    """Select up to a maximum number of the best-seeing images from the specified time range using their Wcs.
    """
    ConfigClass = TimeBestSeeingWcsSelectImageConfig

    def runDataRef(self, dataRef, coordList, makeDataRefList=True,
                   selectDataList=None):
        """Select the best-seeing images in the selectDataList that overlap the patch.
        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Data reference for coadd/tempExp (with tract, patch)
        coordList : `list` of `lsst.afw.geom.SpherePoint`
            List of ICRS sky coordinates specifying boundary of patch
        makeDataRefList : `boolean`, optional
            Construct a list of data references?
        selectDataList : `list` of `SelectStruct`
            List of SelectStruct, to consider for selection
        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:
            - ``exposureList``: the selected exposures
                (`list` of `lsst.pipe.tasks.selectImages.BaseExposureInfo`).
            - ``dataRefList``: the optional data references corresponding to
                each element of ``exposureList``
                (`list` of `lsst.daf.persistence.ButlerDataRef`, or `None`).
        """
        if self.config.nImagesMax and self.config.nImagesMax <= 0:
            raise RuntimeError(f"nImagesMax must be greater than zero: {self.config.nImagesMax}")

        psfSizes = []
        dataRefList = []
        exposureInfoList = []

        if selectDataList is None:
            selectDataList = []

        result = super().runDataRef(dataRef, coordList, makeDataRefList=True, selectDataList=selectDataList)

        for dataRef, exposureInfo in zip(result.dataRefList, result.exposureInfoList):
            cal = dataRef.get("calexp", immediate=True)

            # if min/max PSF values are defined, remove images out of bounds
            psfSize = cal.getPsf().computeShape().getDeterminantRadius()
            sizeFwhm = psfSize * np.sqrt(8.*np.log(2.))

            if self.config.maxPsfFwhm and sizeFwhm > self.config.maxPsfFwhm:
                self.log.debug('Fwhm too large, rejecting, %f', sizeFwhm)
                continue
            if self.config.minPsfFwhm and sizeFwhm < self.config.minPsfFwhm:
                self.log.debug('Fwhm too small, rejecting')
                continue

            mjd = cal.getInfo().getVisitInfo().getDate().get(system=DateTime.MJD)

            if self.config.minMJD and mjd < self.config.minMJD:
                self.log.debug('time too early, rejecting %f', mjd)
                continue

            if self.config.maxMJD and mjd > self.config.maxMJD:
                self.log.debug('time too late, rejecting %f', mjd)
                continue

            psfSizes.append(psfSize)
            dataRefList.append(dataRef)
            exposureInfoList.append(exposureInfo)

        if self.config.nImagesMax and len(psfSizes) > self.config.nImagesMax:
            sortedIndices = np.argsort(psfSizes)[:self.config.nImagesMax]
            filteredDataRefList = [dataRefList[i] for i in sortedIndices]
            filteredExposureInfoList = [exposureInfoList[i] for i in sortedIndices]
            self.log.info(f"{len(sortedIndices)} images selected with FWHM "
                          f"range of {psfSizes[sortedIndices[0]]}--{psfSizes[sortedIndices[-1]]} pixels")

        else:
            if len(psfSizes) == 0:
                self.log.warn(f"0 images selected.")
            else:
                self.log.debug(f"{len(psfSizes)} images selected with FWHM range "
                               f"of {psfSizes[0]}--{psfSizes[-1]} pixels")
            filteredDataRefList = dataRefList
            filteredExposureInfoList = exposureInfoList

        return pipeBase.Struct(
            dataRefList=filteredDataRefList if makeDataRefList else None,
            exposureInfoList=filteredExposureInfoList,
)
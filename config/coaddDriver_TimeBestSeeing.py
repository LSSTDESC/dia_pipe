
from lsst.dia.pipe.selectImages import TimeBestSeeingWcsSelectImagesTask
config.select.retarget(TimeBestSeeingWcsSelectImagesTask)
config.select.minMJD = 59580.0
config.select.maxMJD = 60310.0
config.select.minPsfFwhm = 2.5
config.select.maxPsfFwhm = 3.5
config.select.nImagesMax = 10

from lsst.pipe.drivers.utils import NullSelectImagesTask
config.assembleCoadd.select.retarget(NullSelectImagesTask)

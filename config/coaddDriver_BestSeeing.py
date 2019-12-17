# Load from sub-configurations
import os.path

from lsst.utils import getPackageDir

from lsst.pipe.tasks.assembleCoadd import CompareWarpAssembleCoaddTask
config.assembleCoadd.retarget(CompareWarpAssembleCoaddTask)

for sub, filename in (("makeCoaddTempExp", "makeCoaddTempExp"),
#                      ("backgroundReference", "backgroundReference"),
                      ("assembleCoadd", "compareWarpAssembleCoadd"),
                      ("processCoadd", "processCoadd")):
    path = os.path.join(getPackageDir("obs_lsst"), "config", filename + ".py")
    if os.path.exists(path):
        getattr(config, sub).load(path)

config.doBackgroundReference = False

from lsst.pipe.tasks.selectImages import  BestSeeingWcsSelectImagesTask
config.select.retarget(BestSeeingWcsSelectImagesTask)
config.select.maxPsfFwhm=0.7/0.2

from lsst.pipe.drivers.utils import NullSelectImagesTask
config.assembleCoadd.select.retarget(NullSelectImagesTask)

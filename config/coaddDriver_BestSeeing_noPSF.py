from lsst.pipe.tasks.selectImages import  BestSeeingWcsSelectImagesTask
config.select.retarget(BestSeeingWcsSelectImagesTask)
config.select.maxPsfFwhm=0.7/0.2

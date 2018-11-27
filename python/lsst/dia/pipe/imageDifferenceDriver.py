from lsst.pipe.base import ArgumentParser, ButlerInitializedTaskRunner, ConfigDatasetType
from lsst.pipe.tasks.processCcd import ProcessCcdTask
from lsst.pex.config import Config, Field, ConfigurableField, ListField
from lsst.ctrl.pool.parallel import BatchParallelTask, BatchTaskRunner
from lsst.pipe.tasks.imageDifference import ImageDifferenceTask


class ImageDifferenceDriverConfig(Config):
    ignoreCcdList = ListField(dtype=int, default=[],
                              doc="List of CCDs to ignore when processing")
    ccdKey = Field(dtype=str, default="ccd",
                   doc="DataId key corresponding to a single sensor")
    imageDifference = ConfigurableField(target=ImageDifferenceTask, doc="CCD processing task")



class ImageDifferenceDriverTaskRunner(BatchTaskRunner, ButlerInitializedTaskRunner):
    """Run batches, and initialize Task using a butler"""
    pass


class ImageDifferenceDriverTask(BatchParallelTask):
    """Process CCDs in parallel
    """
    ConfigClass = ImageDifferenceDriverConfig
    _DefaultName = "imageDifferenceDriver"
    RunnerClass = ImageDifferenceDriverTaskRunner

    def __init__(self, butler=None, refObjLoader=None, *args, **kwargs):
        BatchParallelTask.__init__(self, *args, **kwargs)
        self.ignoreCcds = set(self.config.ignoreCcdList)
        self.makeSubtask("imageDifference", butler=butler)

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="imageDifferenceDriver", *args, **kwargs)
        parser.add_id_argument("--id", datasetType="calexp", level="sensor",
                               help="data ID, e.g. --id visit=12345 ccd=67")
        return parser

    def run(self, sensorRef):
        """Process a single CCD, with scatter-gather-scatter using MPI.
        """
        if sensorRef.dataId[self.config.ccdKey] in self.ignoreCcds:
            self.log.warn("Ignoring %s: CCD in ignoreCcdList" %
                          (sensorRef.dataId))
            return None

        with self.logOperation("processing %s" % (sensorRef.dataId,)):
            return self.imageDifference.run(sensorRef)

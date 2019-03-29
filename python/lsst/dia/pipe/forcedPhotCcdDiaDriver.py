from lsst.pipe.base import ArgumentParser, ButlerInitializedTaskRunner
from lsst.pex.config import Config, Field, ConfigurableField, ListField
from lsst.ctrl.pool.parallel import BatchParallelTask, BatchTaskRunner
from .forcedPhotDia import ForcedPhotCcdDiaTask


class ForcedPhotCcdDiaDriverConfig(Config):
    ignoreCcdList = ListField(dtype=int, default=[],
                              doc="List of CCDs to ignore when processing")
    ccdKey = Field(dtype=str, default="detector",
                   doc="DataId key corresponding to a single sensor")
    forced = ConfigurableField(target=ForcedPhotCcdDiaTask, doc="CCD processing task")


class ForcedPhotCcdDiaDriverTaskRunner(BatchTaskRunner, ButlerInitializedTaskRunner):
    """Run batches, and initialize Task using a butler"""
    pass


class ForcedPhotCcdDiaDriverTask(BatchParallelTask):
    """Process CCDs in parallel
    """
    ConfigClass = ForcedPhotCcdDiaDriverConfig
    _DefaultName = "forcedPhotCcdDiaDriver"
    RunnerClass = ForcedPhotCcdDiaDriverTaskRunner

    def __init__(self, butler=None, *args, **kwargs):
        BatchParallelTask.__init__(self, *args, **kwargs)
        self.ignoreCcds = set(self.config.ignoreCcdList)
        self.makeSubtask("forced", butler=butler)

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="forcedPhotCcdDiaDriver", *args, **kwargs)
        parser.add_id_argument("--id", "deepDiff_differenceExp", help="data ID with raw CCD keys"
                               "e.g. --id visit=12345 ccd=1,2")
        return parser

    def runDataRef(self, dataRef):
        """Process a single CCD, with scatter-gather-scatter using MPI.
        """
        if dataRef.dataId[self.config.ccdKey] in self.ignoreCcds:
            self.log.warn("Ignoring %s: CCD in ignoreCcdList" %
                          (dataRef.dataId))
            return None

        with self.logOperation("processing %s" % (dataRef.dataId,)):
            return self.forced.runDataRef(dataRef)

    def writeMetadata(self, dataRef):
        """We don't collect any metadata, so skip"""
        pass

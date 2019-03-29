from lsst.pipe.base import ArgumentParser, TaskRunner
from lsst.pex.config import Config, ConfigurableField, Field
from lsst.pipe.drivers.utils import TractDataIdContainer
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool

from .forcedPhotDia import ForcedPhotCoaddDiaTask


class ForcedPhotCoaddDiaDriverConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")
    forced = ConfigurableField(target=ForcedPhotCoaddDiaTask, doc="Coadd processing task")
    psfCache = Field(dtype=int, default=100, doc="Size of psfCache")


class ForcedPhotCoaddDiaDriverTaskRunner(TaskRunner):
    """Run batches, and initialize Task using a butler"""

    def __init__(self, TaskClass, parsedCmd, doReturnResults=False):
        TaskRunner.__init__(self, TaskClass, parsedCmd, doReturnResults)

    def makeTask(self, parsedCmd=None, args=None):
        """A variant of the base version that passes a butler argument to the task's constructor
        parsedCmd or args must be specified.
        """
        if parsedCmd is not None:
            butler = parsedCmd.butler
        elif args is not None:
            dataRefList, kwargs = args
            butler = dataRefList[0].butlerSubset.butler
        else:
            raise RuntimeError("parsedCmd or args must be specified")
        return self.TaskClass(config=self.config, log=self.log, butler=butler)


def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)


class ForcedPhotCoaddDiaDriverTask(BatchPoolTask):
    """Process Coadds in parallel
    """
    ConfigClass = ForcedPhotCoaddDiaDriverConfig
    _DefaultName = "forcedPhotCoaddDiaDriver"
    RunnerClass = ForcedPhotCoaddDiaDriverTaskRunner

    def __init__(self, butler=None, *args, **kwargs):
        BatchPoolTask.__init__(self, *args, **kwargs)
        self.butler = butler
        self.makeSubtask("forced", butler=butler)

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                                                   parentTask=self._parentTask, log=self.log,
                                                   butler=self.butler))

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="forcedPhotCoaddDiaDriver", *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID with coadd keys"
                               "e.g. --id tract=1234^12345 filter=r^i^z",
                               ContainerClass=TractDataIdContainer)
        return parser

    def runDataRef(self, dataRefList):
        """Process a single Coadd, with scatter-gather-scatter using MPI.
        """
        pool = Pool("all")
        pool.cacheClear()

        pool.map(self.runForced, dataRefList)

    def runForced(self, cache, dataRef):
        """! Run forced photometry on a patch
        Only slave nodes execute this method.
        @param cache: Pool cache, containing butler
        @param patchRef: Patch on which to run measurement
        """
        with self.logOperation("processing %s" % (dataRef.dataId,)):
            return self.forced.runDataRef(dataRef, psfCache=self.config.psfCache)

    def writeMetadata(self, dataRef):
        """We don't collect any metadata, so skip"""
        pass

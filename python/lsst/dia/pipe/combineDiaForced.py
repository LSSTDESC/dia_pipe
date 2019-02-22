import os
import numpy as np
import pandas as pd

import lsst.afw.table as afwTable

from lsst.pex.config import Config, Field, DictField
from lsst.pipe.base import ArgumentParser, TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError
from lsst.pipe.drivers.utils import getDataRef
from lsst.coadd.utils import TractDataIdContainer
from lsst.daf.base import DateTime


__all__ = ("CombineDiaForcedConfig", "CombineDiaForcedTask")


class CombineDiaForcedConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")
    ccdKey = Field(dtype=str, default='detector', doc="Name of ccd to give to butler")
    keepFields = DictField(
        keytype=str, itemtype=str,
        doc='Keep these fields from the diaSrc catalogs.  Will be averaged over matches',
        default={
            'psf_flux': 'base_PsfFlux_instFlux',
            'psf_flux_err': 'base_PsfFlux_instFluxErr'
        }
    )
    storage = Field(dtype=str, default="pickle", doc="pandas storage format")


class CombineDiaForcedTaskRunner(TaskRunner):
    """TaskRunner for running CombineDiaForced Task.

    This is similar to the lsst.pipe.drivers.multiBandDriver runner
    except with no reuse outputs option.
    """

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

class CombineDiaForcedTask(BatchPoolTask):
    """Multi-node driver for multiband processing"""
    ConfigClass = CombineDiaForcedConfig
    _DefaultName = "CombineDiaForced"
    RunnerClass = CombineDiaForcedTaskRunner

    def __init__(self, butler=None, **kwargs):
        """!
        @param[in] butler: the butler can be used to retrieve schema or passed to the refObjLoader
            constructor in case it is needed.
        """
        BatchPoolTask.__init__(self, **kwargs)
        self.butler = butler
        self.schema = afwTable.SourceTable.makeMinimalSchema()

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                                                   parentTask=self._parentTask, log=self.log,
                                                   butler=self.butler))

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name=cls._DefaultName, *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd",
                               help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=TractDataIdContainer)
        return parser

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numCpus):
        """!Return walltime request for batch job

        @param time: Requested time per iteration
        @param parsedCmd: Results of argument parsing
        @param numCores: Number of cores
        """
        numTargets = 0
        for refList in parsedCmd.id.refList:
            numTargets += len(refList)
        return time*numTargets/float(numCpus)

    @abortOnError
    def runDataRef(self, patchRefList):
        """!Combine forced diaObjects into a single catalog to construct light curves

        Only the master node runs this method.

        @param patchRefList:  Data references to run measurement
        """
        for patchRef in patchRefList:
            if patchRef:
                butler = patchRef.getButler()
                break
        else:
            raise RuntimeError("No valid patches")
        pool = Pool("all")
        pool.cacheClear()
        pool.storeSet(butler=butler)

        # Group all filters by patch
        patches = {}
        tract = None
        for patchRef in patchRefList:
            dataId = patchRef.dataId
            if tract is None:
                tract = dataId["tract"]
            else:
                assert tract == dataId["tract"]

            patch = dataId["patch"]
            if patch not in patches:
                patches[patch] = []
            patches[patch].append(dataId)

        pool.map(self.runCombine, patches.values())

    def runCombine(self, cache, dataIdList):
        """! Run combination on a patch

        For all of the visits that overlap this patch in the band create a catalog with
        mjd, id, filter and values in the keepFields config parameter.

        It will write a file as a pandas DataFrame to the appropriate deepDiff directory,
        but the products are not currently declared in the mapper file.
        """

        dataRefList = [getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd_calexp") for
                       dataId in dataIdList]

        try:
            diaObject = dataRefList[0].get(f"{self.config.coaddName}Diff_diaObject")
            uri = dataRefList[0].getUri(f"{self.config.coaddName}Diff_diaObject")
        except Exception:
            self.log.info('Cannot read diaObject for %s' % (dataRefList[0].dataId))
            return

        # Use a dictionary of arrays to store data.  Probably something more efficient
        data = {}
        data['id'] = []
        data['mjd'] = []
        data['filter'] = []
        for key, value in self.config.keepFields.items():
            data[key] = []

        for dataRef in dataRefList:

            try:
                calexp = dataRef.get(f"{self.config.coaddName}Coadd_calexp")
            except Exception:
                self.log.info('Cannot read data for %s' % (dataRef.dataId))
                continue

            visitCatalog = calexp.getInfo().getCoaddInputs().ccds

            for visitRec in visitCatalog:

                visit = int(visitRec.get('visit'))
                ccd = int(visitRec.get('ccd'))
                dataId = {"visit": visit, self.config.ccdKey: ccd}

                try:
                    src = cache.butler.get(f"{self.config.coaddName}Diff_forced_dia_src", dataId)
                    diff = cache.butler.get(f"{self.config.coaddName}Diff_differenceExp", dataId)
                except Exception as e:
                        self.log.debug('Cannot read data for %d %d. skipping %s', visit, ccd, e)
                        continue

                mjd = diff.getInfo().getVisitInfo().getDate().get(system=DateTime.MJD)
                band = diff.getInfo().getFilter().getName()
                self.log.info('Reading diff forced src with %d sources %s', len(src), dataId)

                matches = np.in1d(src['dia_object_id'], diaObject['id'], assume_unique=True)

                data['id'].extend(src['dia_object_id'][matches])
                data['mjd'].extend([mjd]*np.sum(matches))
                data['filter'].extend([band]*np.sum(matches))
                for key, value in self.config.keepFields.items():
                    data[key].extend(src[value][matches])

        tract = dataRefList[0].dataId['tract']
        patch = dataRefList[0].dataId['patch']
        path = os.path.dirname(uri)
        df = pd.DataFrame(data)
        getattr(df, "to_" + self.config.storage)(f'{path}/diaCombined_{tract}_{patch}.{self.config.storage}')

    def writeMetadata(self, dataRef):
        """We don't collect any metadata, so skip"""
        pass

    def _getConfigName(self):
        """!Return the name of the config dataset.
        """
        pass

    def _getMetadataName(self):
        """!Return the name of the metadata dataset.
        """
        pass

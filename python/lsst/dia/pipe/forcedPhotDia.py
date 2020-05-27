import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig

from lsst.meas.base.references import MultiBandReferencesConfig, BaseReferencesTask, CoaddSrcReferencesTask
from lsst.meas.base.forcedMeasurement import ForcedPluginConfig, ForcedPlugin
from lsst.meas.base.pluginRegistry import register


class DiaReferencesConfig(CoaddSrcReferencesTask.ConfigClass):
    datasetSuffix = pexConfig.Field(dtype=str, default='deepDiff_diaObject',
                                    doc="Name of dataset to get")
    def validate(self):
        if self.filter is not None:
            raise pexConfig.FieldValidationError(
                field=MultiBandReferencesConfig.filter,
                config=self,
                msg="Filter should not be set for the multiband processing scheme")
        # Delegate to ultimate base class, because the direct one has a check we don't want.
        BaseReferencesTask.ConfigClass.validate(self)


class DiaReferencesTask(CoaddSrcReferencesTask):
    """!
    Dummy task so that the reference schema will match

    We do not get the reference tracts/patches using this class
    """
    ConfigClass = DiaReferencesConfig

    def __init__(self, butler=None, schema=None, **kwargs):
        """! Initialize the task.
        We cannot use the CoaddSrcReferencesTask init because of the hard coded "Coadd"

        Additional keyword arguments (forwarded to BaseReferencesTask.__init__):
         - schema: the schema of the detection catalogs used as input to this one
         - butler: a butler used to read the input schema from disk, if schema is None
        The task will set its own self.schema attribute to the schema of the output merged catalog.
        """
        BaseReferencesTask.__init__(self, butler=butler, schema=schema, **kwargs)
        if schema is None:
            assert butler is not None, "No butler nor schema provided"
            schema = butler.get("{}_schema".format(self.config.datasetSuffix),
                                immediate=True).getSchema()
        self.schema = schema



class DiaSrcReferencesConfig(BaseReferencesTask.ConfigClass):
    datasetSuffix = pexConfig.Field(dtype=str, default='deepDiff_diaSrc',
                                    doc="Name of dataset to get")


class DiaSrcReferencesTask(BaseReferencesTask):
    """!
    Dummy task so that the reference schema will match

    We do not get the reference tracts/patches using this class
    """
    ConfigClass = DiaSrcReferencesConfig

    def __init__(self, butler=None, schema=None, **kwargs):
        """! Initialize the task.
        We cannot use the CoaddSrcReferencesTask init because of the hard coded "Coadd"

        Additional keyword arguments (forwarded to BaseReferencesTask.__init__):
         - schema: the schema of the detection catalogs used as input to this one
         - butler: a butler used to read the input schema from disk, if schema is None
        The task will set its own self.schema attribute to the schema of the output merged catalog.
        """
        BaseReferencesTask.__init__(self, butler=butler, schema=schema, **kwargs)
        if schema is None:
            assert butler is not None, "No butler nor schema provided"
            schema = butler.get("{}_schema".format(self.config.datasetSuffix),
                                immediate=True).getSchema()
        self.schema = schema



class ForcedDiaTransformedCentroidConfig(ForcedPluginConfig):
    pass


@register("base_DiaTransformedCentroid")
class ForcedDiaTransformedCentroidPlugin(ForcedPlugin):
    """A centroid pseudo-algorithm for forced measurement that simply transforms sky coordinate
    from the reference catalog to the measurement coordinate system.  This does assumes that
    there is no centroid column in the refernce catalog.
    """

    ConfigClass = ForcedDiaTransformedCentroidConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        # Allocate x and y fields, join these into a single FunctorKey for ease-of-use.
        xKey = schema.addField(name + "_x", type="D", doc="transformed reference centroid column",
                               units="pixel")
        yKey = schema.addField(name + "_y", type="D", doc="transformed reference centroid row",
                               units="pixel")
        self.centroidKey = afwTable.Point2DKey(xKey, yKey)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        targetWcs = exposure.getWcs()
        targetPos = targetWcs.skyToPixel(refRecord.getCoord())
        measRecord.set(self.centroidKey, targetPos)


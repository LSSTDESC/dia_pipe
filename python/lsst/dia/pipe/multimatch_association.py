"""MultiMatch sequenctial association of DIASources into DIAObjects.
"""
import numpy as np
from collections import defaultdict
import pandas as pd

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom as geom

from .parquetTable import ParquetTable

__all__ = ["MultiMatchAssociationConfig", "MultiMatchAssociationTask"]


class MultiMatchAssociationConfig(pexConfig.Config):
    """Configuration parameters for the MultiMatchAssociationTask
    """
    tolerance = pexConfig.Field(
        dtype=float,
        doc='maximum distance to match sources together in arcsec',
        default=0.5
    )
    fluxType = pexConfig.Field(
        dtype=str,
        doc='Keep track of the average flux of this type',
        default='base_PsfFlux_instFlux',
    )


class MultiMatchAssociationTask(pipeBase.Task):
    """Construct DIAObjects from a list of DIASources
    """

    ConfigClass = MultiMatchAssociationConfig
    _DefaultName = "MultiMatch_association"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.multi_matches = None

    def addCatalog(self, src, filt, visit, ccd, footprints):
        """Add objects from a catalog to the existing MultiMatch

        @param[in]  srcCat      An SourceCatalog of objects to be added.
        @param[in]  filt        The filter of the catalog
        @param[in]  visit       The visit number
        @param[in]  ccd         The ccd number
        @param[in]  footprints  A list of footprints that have been transformed to the
                                WCS of the coadd patch.
        """

        if self.multi_matches is None:
            self.multi_matches = afwTable.MultiMatch(src.schema, {'visit': np.int32, 'ccd': np.int32},
                                                     radius=afwGeom.Angle(self.config.tolerance/3600.,
                                                     geom.degrees))
        for s, foot in zip(src, footprints):
            s.setFootprint(foot)
        self.multi_matches.add(src, {'visit': visit, 'ccd': ccd})

    def initialize(self, schema, idFactory):
        pass

    def finalize(self, idFactory):
        """Finalize construction of the catalog.

        Create a SourceCatalog from the MultiMatch object and compute the corresponding
        merged footprint.
        @param[in]  idFactory   Used to generate ids.

        @return    SourceCatalog of DIAObjects
        """
        if self.multi_matches is None:
            return None

        schema = afwTable.SourceTable.makeMinimalSchema()
        nobsKey = schema.addField("nobs", type=np.int32, doc='Number of times observed')
        fluxKey = schema.addField("flux", type=float, doc='Average flux')
        raKey = schema['coord_ra'].asKey()
        decKey = schema['coord_dec'].asKey()
        table = afwTable.SourceTable.make(schema, idFactory)
        cat = afwTable.SourceCatalog(table)

        results = self.multi_matches.finish(removeAmbiguous=False)
        allMatches = afwTable.GroupView.build(results)

        psfMagKey = allMatches.schema.find(self.config.fluxType).key
        raKey = allMatches.schema.find("coord_ra").key
        decKey = allMatches.schema.find("coord_dec").key

        ave_ra = allMatches.aggregate(np.mean, field=raKey)
        ave_dec = allMatches.aggregate(np.mean, field=decKey)
        ave_flux = allMatches.aggregate(np.mean, field=psfMagKey)

        # Merge the footprints from the same object together
        object_ids = np.unique(results['object'])
        footprints = []
        self.diaSrcIds = []
        self.diaObjectIds = []
        for id in object_ids:
            mask = results['object'] == id
            footprint = None
            src_ids = []
            for rec in results[mask]:
                if footprint is None:
                    footprint = rec.getFootprint()
                else:
                    footprint = afwDet.mergeFootprints(footprint, rec.getFootprint())
                src_ids.append(rec.get('id'))
            self.diaSrcIds.append(src_ids)
            footprints.append(footprint)

        for i in range(len(ave_ra)):
            rec = cat.addNew()
            self.diaObjectIds.append(rec.get('id'))
            rec.setFootprint(footprints[i])
            rec.set(raKey, ave_ra[i]*geom.radians)
            rec.set(decKey, ave_dec[i]*geom.radians)
            rec.set(fluxKey, ave_flux[i])

        return cat

    def getObjectIds(self):
        """Get a list of id's corresponding to the objects in this catalog

        @return    pandas DataFrame of matching diaObject ids to diaSrc ids
        """
        data = {}
        data['diaObjectId'] = self.diaObjectIds
        data['diaSrcIds'] = self.diaSrcIds
        df = pd.DataFrame(data)
        table = ParquetTable(dataFrame=df)
        return table

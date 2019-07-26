"""Simple sequenctial association of DIASources into DIAObjects.
"""
import numpy as np
import pandas as pd

import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from .utils import query_disc, eq2xyz, toIndex
from .parquetTable import ParquetTable

__all__ = ["SimpleAssociationConfig", "SimpleAssociationTask"]


class SimpleAssociationConfig(pexConfig.Config):
    """Configuration parameters for the SimpleAssociationTask
    """
    tolerance = pexConfig.Field(
        dtype=float,
        doc='maximum distance to match sources together in arcsec',
        default=0.5
    )
    nside = pexConfig.Field(
        dtype=int,
        doc='Healpix nside value used for indexing',
        default=2**18,
    )
    aveFields = pexConfig.ListField(
        dtype=str,
        doc='Average these fields from the diaSrc catalogs per band.  '
            'They must have a corresponding Err quantity.',
        default=['base_PsfFlux_instFlux']
    )
    filters = pexConfig.ListField(
        dtype=str,
        doc='Which filters will be averaged over',
        default=['u', 'g', 'r', 'i', 'z', 'y']
    )


class SimpleAssociationTask(pipeBase.Task):
    """Construct DIAObjects from a list of DIASources
    """

    ConfigClass = SimpleAssociationConfig
    _DefaultName = "simple_association"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.cat = None
        self.footprints = []

    def dist(self, src_ra, src_dec, tol):
        """Compute the distance to all the objects within the tolerance of a point.

        @param[in]  src_ra      RA in degrees
        @param[in]  src_dec     Dec in degrees
        @param[in]  tol         Matching tolerance in arcseconds

        @return    Pair of distances and matching indices to existing catalog
        """
        ra = np.rad2deg(self.cat[self.keys['coord_ra']])
        dec = np.rad2deg(self.cat[self.keys['coord_dec']])
        match_indices = query_disc(self.config.nside, src_ra, src_dec, np.deg2rad(tol/3600.))
        matches = np.in1d(self.indexes, match_indices)

        if np.sum(matches) < 1:
            return None, None

        dist = np.array([np.sqrt(np.sum((eq2xyz(src_ra, src_dec) - eq2xyz(ra1, dec1))**2))
                         for ra1, dec1 in zip(ra[matches], dec[matches])])
        return dist, matches

    def initialize(self, src_schema, idFactory):
        """Initialize the catalog

        We need to know the schema before constructing the catalog, and thus we also need
        an id factory.

        @param[in]  src_schema  Schema from the diaSource
        @param[in]  idFactory   Used to generate ids.
        """
        self.schema = afwTable.SourceTable.makeMinimalSchema()

        self.keys = {}
        self.keys['nobs'] = self.schema.addField("nobs", type=np.int32,
                                                 doc='Number of times observed in all filters')
        self.keys['coord_ra'] = self.schema['coord_ra'].asKey()
        self.keys['coord_dec'] = self.schema['coord_dec'].asKey()

        new_fields = {"Mean", "MeanErr", "Sigma", "Chi2", "Ndata"}

        for filter in self.config.filters:
            for field_name in self.config.aveFields:
                for new_field in new_fields:
                    field = src_schema[field_name].asField()
                    name = f"{field.getName()}_{new_field}_{filter}"

                    if new_field == "Ndata":
                        field_type = np.int32
                        doc = f'Number of times observed in filter {filter}'
                    else:
                        field_type = field.getTypeString()
                        doc = f"{new_field} {field.getDoc()} for the {filter} filter"
                    added_field = self.schema.addField(name, type=field_type, doc=doc)
                    self.keys[name] = added_field

        self.table = afwTable.SourceTable.make(self.schema, idFactory)
        self.cat = afwTable.SourceCatalog(self.table)

        # store the healpix index
        self.indexes = []
        self.idLists = []
        self.footprints = []

    def addNew(self, src, footprint, filter):
        """Add a new object to the existing catalog.

        @param[in]  src         SourceRecofd of new object
        @param[in]  footprint   The object footprint in the coadd WCS patch.
        @param[in]  filter   The object filter
        """
        rec = self.cat.addNew()
        rec.set(self.keys['nobs'], 1)
        rec.set(self.keys['coord_ra'], src.get('coord_ra'))
        rec.set(self.keys['coord_dec'], src.get('coord_dec'))

        for field in self.config.aveFields:
            val = src.get(field)
            if np.isfinite(val):
                rec.set(self.keys[f"{field}_Mean_{filter}"], val)
                rec.set(self.keys[f"{field}_Sigma_{filter}"], 0.)
                rec.set(self.keys[f"{field}_MeanErr_{filter}"], 0.)
                rec.set(self.keys[f"{field}_Chi2_{filter}"], 0.)
                rec.set(self.keys[f"{field}_Ndata_{filter}"], 1)

        index = toIndex(self.config.nside, src.get('coord_ra').asDegrees(),
                        src.get('coord_dec').asDegrees())
        self.indexes.append(index)
        self.idLists.append([src.get('id')])
        self.footprints.append(footprint)

    def updateCat(self, match, src, footprint, filter):
        """Update an object in the existing catalog.

        @param[in]  match       Index of the current catalog that matches
        @param[in]  src         SourceRecord of the matching object
        @param[in]  footprint   The object footprint in the coadd WCS patch.
        @param[in]  filter      The filter of the catalog
        """
        self.idLists[match].append(src.get('id'))
        match = int(match)

        nobs = self.cat[match].get(self.keys['nobs'])
        self.cat[match].set(self.keys['nobs'], nobs + 1)

        # update coordinates as average over all bands
        ra = self.cat[match].get(self.keys['coord_ra'])
        dec = self.cat[match].get(self.keys['coord_dec'])

        self.cat[match].set(self.keys['coord_ra'],
                            (src.get('coord_ra') + (nobs - 1)*ra)/nobs)
        self.cat[match].set(self.keys['coord_dec'],
                            (src.get('coord_dec') + (nobs - 1)*dec)/nobs)

        self.footprints[match] = afwDet.mergeFootprints(self.footprints[match], footprint)

        for field in self.config.aveFields:
            new_val = src.get(field)

            if np.isfinite(new_val) is False:
                continue

            ndata = self.cat[match].get(f"{field}_Ndata_{filter}")
            # If object was first detected in another filter then it can still
            # have zero entries for other filters
            if ndata == 0:
                self.cat[match].set(self.keys[f"{field}_Mean_{filter}"], new_val)
                self.cat[match].set(self.keys[f"{field}_Sigma_{filter}"], 0.)
                self.cat[match].set(self.keys[f"{field}_MeanErr_{filter}"], 0.)
                self.cat[match].set(self.keys[f"{field}_Chi2_{filter}"], 0.)
                self.cat[match].set(self.keys[f"{field}_Ndata_{filter}"], 1)
            else:
                ndata += 1
                mean_val = self.cat[match].get(f"{field}_Mean_{filter}")
                sigma_val = self.cat[match].get(f"{field}_Sigma_{filter}")
                chi2_val = self.cat[match].get(f"{field}_Chi2_{filter}")

                new_val_error = src.get(field + "Err")

                new_mean = (new_val + (ndata - 1)*mean_val)/ndata
                new_sigma = np.sqrt(sigma_val +
                                    (new_val - mean_val)*(new_val - new_mean))/ndata
                new_chi2 = chi2_val + ((new_val - mean_val)/new_val_error)**2

                self.cat[match].set(self.keys[f"{field}_Mean_{filter}"], new_mean)
                self.cat[match].set(self.keys[f"{field}_Sigma_{filter}"], new_sigma)
                self.cat[match].set(self.keys[f"{field}_MeanErr_{filter}"], new_sigma/np.sqrt(ndata))
                self.cat[match].set(self.keys[f"{field}_Chi2_{filter}"], new_chi2)
                self.cat[match].set(self.keys[f"{field}_Ndata_{filter}"], ndata)

    def addCatalog(self, srcCat, filter, visit, ccd, footprints):
        """Add objects from a new catalog to the existing list

        For objects that are not within the tolerance of any existing objects, new
        sources will be created, otherwise the existing objects will be update.
        @param[in]  srcCat      An SourceCatalog of objects to be added.
        @param[in]  filt        The filter of the catalog
        @param[in]  visit       The visit number
        @param[in]  ccd         The ccd number
        @param[in]  footprints  A list of footprints that have been transformed to the
                                WCS of the coadd patch.
        """
        for src, footprint in zip(srcCat, footprints):
            dist, matches = self.dist(src['coord_ra'].asDegrees(), src['coord_dec'].asDegrees(),
                                      2*self.config.tolerance)

            if dist is None:
                self.addNew(src, footprint, filter)
                continue

            if np.min(dist) < np.deg2rad(self.config.tolerance/3600):
                match_dist = np.argmin(dist)
                match_index = np.where(matches)[0][match_dist]
                self.updateCat(match_index, src, footprint, filter)
            else:
                self.addNew(src, footprint, filter)

    def finalize(self, idFactory):
        """Finalize construction by setting the footprint

        @param[in]  idFactory   Used to generate ids.

        @return    SourceCatalog of DIAObjects
        """
        if self.cat is None:
            return None

        for rec, footprint in zip(self.cat, self.footprints):
            rec.setFootprint(footprint)

        return self.cat

    def getObjectIds(self):
        """Get a list of id's corresponding to the objects in this catalog

        @return    pandas DataFrame of matching diaObject ids to diaSrc ids
        """
        data = {}
        data['diaObjectId'] = self.cat['id']
        data['diaSrcIds'] = self.idLists
        df = pd.DataFrame(data)
        table = ParquetTable(dataFrame=df)
        return table


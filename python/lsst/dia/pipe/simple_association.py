"""Simple sequenctial association of DIASources into DIAObjects.
"""
import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from .utils import query_disc, eq2xyz, toIndex

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
    keepFields = pexConfig.ListField(
        dtype=str,
        doc='Keep these fields from the diaSrc catalogs.  Will be averaged over matches',
        default=['base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']
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
        self.keys['nobs'] = self.schema.addField("nobs", type=np.int32, doc='Number of times observed')
        self.keys['coord_ra'] = self.schema['coord_ra'].asKey()
        self.keys['coord_dec'] = self.schema['coord_dec'].asKey()

        for field in self.config.keepFields:
            self.keys[field] = self.schema.addField(src_schema[field].asField())

        self.table = afwTable.SourceTable.make(self.schema, idFactory)
        self.cat = afwTable.SourceCatalog(self.table)
        # store the healpix index
        self.indexes = []
        self.idLists = []
        self.footprints = []

    def addNew(self, src, footprint):
        """Add a new object to the existing catalog.

        @param[in]  src         SourceRecofd of new object
        @param[in]  footprint   The object footprint in the coadd WCS patch.
        """
        rec = self.cat.addNew()
        rec.set(self.keys['nobs'], 1)
        for label, key in self.keys.items():
            if label in src.schema:
                rec.set(key, src.get(label))
        index = toIndex(self.config.nside, src.get('coord_ra').asDegrees(),
                        src.get('coord_dec').asDegrees())
        self.indexes.append(index)
        self.idLists.append([src.get('id')])
        self.footprints.append(footprint)

    def updateCat(self, match, src, footprint):
        """Update an object in the existing catalog.

        @param[in]  match       Index of the current catalog that matches
        @param[in]  src         SourceRecord of the matching object
        @param[in]  footprint   The object footprint in the coadd WCS patch.
        """
        self.idLists[match].append(src.get('id'))
        match = int(match)
        nobs = self.cat[match].get(self.keys['nobs'])
        self.cat[match].set(self.keys['nobs'], nobs + 1)

        self.footprints[match] = afwDet.mergeFootprints(self.footprints[match], footprint)

        for label, key in self.keys.items():
            if label in src.schema:
                val = self.cat[match].get(label)
                new_val = src.get(label)
                self.cat[match].set(key, ((nobs - 1)*val + new_val)/nobs)

    def addCatalog(self, srcCat, filt, visit, ccd, footprints):
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
                self.addNew(src, footprint)
                continue

            if np.min(dist) < np.deg2rad(self.config.tolerance/3600):
                match_dist = np.argmin(dist)
                match_index = np.where(matches)[0][match_dist]
                self.updateCat(match_index, src, footprint)
            else:
                self.addNew(src, footprint)

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

"""Simple sequenctial association of DIASources into DIAObjects.
"""

import numpy as np
import healpy as hp

import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom as geom

from .utils import query_disc, eq2xyz, toIndex

__all__ = ["SimpleAssociationConfig", "SimpleAssociationTask"]

class DIAObjectList:

	def __init__(self, nside, ra=None, dec=None, idLists=None, nobs=None, indexes=None, 
				 filters=None, fluxes=None, footprints=None):
		self.nside = nside
		if ra is not None:
			if np.all([ra.size == dec.size, ra.size == len(idLists), 
					   ra.size == nobs.size, ra.size == indexes.size,
					   ra.size == len(filters), ra.size == fluxes.size,
					   ra.size == len(footprints)]) is False:
				raise ValueError("DIAObjectList arrays must all be the same size")
			self.ra = ra
			self.dec = dec
			self.idLists = idLists
			self.nobs = nobs
			self.indexes = indexes
			self.filters = filters
			self.fluxes = fluxes
			self.footprints = footprints
		else:
			self.ra = np.array([])
			self.dec = np.array([])
			self.idLists = []            
			self.nobs = np.array([])
			self.filters = []
			self.fluxes = np.array([])
			self.indexes = np.array([])
			self.footprints = []

	def dist(self, ra, dec, tol):
		match_indices = query_disc(self.nside, ra, dec, np.deg2rad(tol/3600))
		matches = np.in1d(self.indexes, match_indices)
		if np.sum(matches) < 1:
		    return None, None
		dist = np.array([np.sqrt(np.sum((eq2xyz(ra,dec) - eq2xyz(ra1,dec1))**2)) 
						 for ra1,dec1 in zip(self.ra[matches], self.dec[matches])])
		return dist, matches

	def update(self, match, ra, dec, id, filt, flux, footprint):
		self.idLists[match].append(id)
		self.filters[match].append(filt)
		self.nobs[match] += 1
		nobs = self.nobs[match]
		self.ra[match] = ((nobs - 1)*self.ra[match] + ra)/nobs
		self.dec[match] = ((nobs - 1)*self.dec[match] + dec)/nobs
		self.fluxes[match] = ((nobs - 1)*self.fluxes[match] + flux)/nobs
		self.footprint[match] = afwDet.mergeFootprints(self.footprint[match], footprint)

	def extend(self, other):
		self.indexes = np.append(self.indexes, other.indexes)
		self.ra = np.append(self.ra, other.ra)
		self.dec = np.append(self.dec, other.dec)
		self.nobs = np.append(self.nobs, other.nobs)
		self.idLists.extend(other.idLists)
		self.filters.extend(other.filters)
		self.fluxes = np.append(self.fluxes, other.fluxes)
		self.footprints.extend(other.footprints)

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
    fluxType = pexConfig.Field(
        dtype=str,
        doc='Keep track of the average flux of this type',
        default='base_PsfFlux_flux',
    )


class SimpleAssociationTask(pipeBase.Task):
	"""Construct DIAObjects from a list of DIASources
	"""

	ConfigClass = SimpleAssociationConfig
	_DefaultName = "simple_association"

	def __init__(self, **kwargs):

		pipeBase.Task.__init__(self, **kwargs)
		self.dia_cat = DIAObjectList(self.config.nside)
		self.cat = None
		

	def addCatalog(self, src, filt, visit, ccd, footprints):

		src_ra = np.rad2deg(src['coord_ra'])
		src_dec = np.rad2deg(src['coord_dec'])
		src_id = src['id']
		src_fluxes = src[self.config.fluxType]
		src_indexes = toIndex(self.config.nside, src_ra, src_dec)
		new_ra = []
		new_dec = []
		new_index = []
		new_ids = []
		new_nobs = []
		new_filters = []
		new_fluxes = []
		new_footprints = []

		for ii,(ra,dec,index,id,flux,footprint) in enumerate(zip(src_ra, src_dec, 
																 src_indexes, src_id, 
																 src_fluxes, footprints)):
			dist, matches = self.dia_cat.dist(ra, dec, 2*self.config.tolerance)
			# Create a new object to be added 
			if dist is None:
			    new_ra.append(ra)
			    new_dec.append(dec)
			    new_index.append(index)
			    new_ids.append([id])
			    new_nobs.append(1)
			    new_filters.append(filt)
			    new_fluxes.append(flux)
			    new_footprints.append(footprint)
			continue

			if np.min(dist) < np.deg2rad(self.config.tolerance/3600):
			    match_dist = np.argmin(dist)
			    match_index = np.where(matches)[0][match_dist]
			    self.dia_cat.update(match_index, ra, dec, index, filt, flux, footprint)
			else:
			    new_ra.append(ra)
			    new_dec.append(dec)
			    new_index.append(index)
			    new_ids.append([id])
			    new_nobs.append(1)
			    new_filters.append(filt)
			    new_fluxes.append(flux)
			    new_footprints.append(footprint)

		if len(new_ra) > 0:
			new_cat = DIAObjectList(self.config.nside, ra=np.array(new_ra), dec=np.array(new_dec), 
									indexes=np.array(new_index), nobs=np.array(new_nobs), 
									idLists=new_ids, filters=new_filters,
									fluxes=np.array(new_fluxes), footprints=new_footprints)
			self.log.info('Adding %d new objects' % len(new_ra))
			self.dia_cat.extend(new_cat)

	def finalize(self, idFactory):
		"""Finalize construction by creating afwTable SourceCatalog"""
		schema = afwTable.SourceTable.makeMinimalSchema()
		nobsKey = schema.addField("nobs", type=np.int32, doc='Number of times observed')
		fluxKey = schema.addField("flux", type=float, doc='Average flux')
		raKey = schema['coord_ra'].asKey()
		decKey = schema['coord_dec'].asKey()
		table = afwTable.SourceTable.make(schema, idFactory)
		cat = afwTable.SourceCatalog(table)

		for i in range(len(self.dia_cat.ra)):
			rec = cat.addNew()
			rec.setFootprint(self.dia_cat.footprints[i])
			rec.set(raKey, self.dia_cat.ra[i]*geom.degrees)
			rec.set(decKey, self.dia_cat.dec[i]*geom.degrees)
			rec.set(fluxKey, self.dia_cat.fluxes[i])
			rec.set(nobsKey, int(self.dia_cat.nobs[i]))

		return cat



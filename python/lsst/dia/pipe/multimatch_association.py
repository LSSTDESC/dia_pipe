"""MultiMatch sequenctial association of DIASources into DIAObjects.
"""
import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom as geom

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
        default='base_PsfFlux_flux',
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

		if self.multi_matches is None:
			self.multi_matches = afwTable.MultiMatch(src.schema, {'visit':np.int32, 'ccd':np.int32}, 
				                					 radius=afwGeom.Angle(self.config.tolerance/3600., geom.degrees))
		for s,foot in zip(src, footprints):
			s.setFootprint(foot)
		self.multi_matches.add(src, {'visit':visit, 'ccd':ccd})

	def finalize(self, idFactory):
		"""Finalize construction by creating afwTable SourceCatalog"""

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
		for id in object_ids:
			mask = results['object'] == id
			footprint = None
			for rec in results[mask]:
				if footprint is None:
					footprint = rec.getFootprint()
				else:
					footprint = afwDet.mergeFootprints(footprint, rec.getFootprint())
			footprints.append(footprint)


		for i in range(len(ave_ra)):
			rec = cat.addNew()
			rec.setFootprint(footprints[i])
			rec.set(raKey, ave_ra[i]*geom.radians)
			rec.set(decKey, ave_dec[i]*geom.radians)
			rec.set(fluxKey, ave_flux[i])

		return cat



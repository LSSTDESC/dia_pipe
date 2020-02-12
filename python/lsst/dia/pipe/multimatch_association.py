"""MultiMatch sequenctial association of DIASources into DIAObjects.
"""
import numpy as np
from collections import defaultdict
import pandas as pd

import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom as geom
import lsst.afw.image as afwImage

from .parquetTable import ParquetTable

__all__ = ["MultiMatchAssociationConfig", "MultiMatchAssociationTask"]


def scaleFlux(flux, flux_err, calib, new_calib):
    """Scale flux and error to new zeropoint
    """
    mag = calib.instFluxToMagnitude(flux, flux_err)
    flux = new_calib.magnitudeToInstFlux(mag.value)
    flux_err = flux*0.4*np.log(10)*mag.error
    return flux, flux_err


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
    filters = pexConfig.ListField(
        dtype=str,
        doc='Which filters will be averaged over',
        default=['u', 'g', 'r', 'i', 'z', 'y']
    )
    commonZp = pexConfig.Field(
        dtype=float,
        doc='Put all fluxes on common zeropoint',
        default=27
    )


class MultiMatchAssociationTask(pipeBase.Task):
    """Construct DIAObjects from a list of DIASources
    """

    ConfigClass = MultiMatchAssociationConfig
    _DefaultName = "MultiMatch_association"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.multi_matches = None
        self.calib = afwImage.makePhotoCalibFromCalibZeroPoint(10**(0.4*self.config.commonZp))
        self.calibDict = defaultdict(dict)

    def addCatalog(self, src, filter, visit, ccd, calib, footprints):
        """Add objects from a catalog to the existing MultiMatch

        @param[in]  srcCat      An SourceCatalog of objects to be added.
        @param[in]  filt        The filter of the catalog
        @param[in]  visit       The visit number
        @param[in]  ccd         The ccd number
        @param[in]  footprints  A list of footprints that have been transformed to the
                                WCS of the coadd patch.
        """

        if self.multi_matches is None:
            # The data id for multiMatch does not take strings so we need to convert filter to a string
            self.multi_matches = afwTable.MultiMatch(src.schema, {'visit': np.int32, 'ccd': np.int32,
                                                                  'filter': np.int32},
                                                     radius=geom.Angle(self.config.tolerance/3600.,
                                                     geom.degrees))
        for s, foot in zip(src, footprints):
            s.setFootprint(foot)
        self.multi_matches.add(src, {'visit': visit, 'ccd': ccd,
                                     'filter': self.config.filters.index(filter)})
        self.calibDict[visit][ccd] = calib

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

        keys = {}
        for filter in self.config.filters:
            flux = self.config.fluxType
            keys[f'{flux}_Mean_{filter}'] = schema.addField(f"{flux}_Mean_{filter}", type=float,
                                                            doc=f'Mean {flux} in filter {filter}')
            keys[f'{flux}_MeanErr_{filter}'] = schema.addField(f"{flux}_MeanErr_{filter}", type=float,
                                                               doc=f'MeanErr {flux} in filter {filter}')
            keys[f'{flux}_Sigma_{filter}'] = schema.addField(f"{flux}_Sigma_{filter}", type=float,
                                                             doc=f'Sigma {flux} in filter {filter}')
            keys[f'{flux}_Ndata_{filter}'] = schema.addField(f"{flux}_NData_{filter}", type=np.int32,
                                                             doc=f'Number of observations in filter {filter}')
            keys[f'{flux}_Chi2_{filter}'] = schema.addField(f"{flux}_Chi2_{filter}", type=float,
                                                            doc=f'Chi2 of {flux} for {filter}')
        raKey = schema['coord_ra'].asKey()
        decKey = schema['coord_dec'].asKey()
        table = afwTable.SourceTable.make(schema, idFactory)
        cat = afwTable.SourceCatalog(table)

        results = self.multi_matches.finish(removeAmbiguous=False)
        allMatches = afwTable.GroupView.build(results)

        raKey = allMatches.schema.find("coord_ra").key
        decKey = allMatches.schema.find("coord_dec").key

        ave_ra = allMatches.aggregate(np.mean, field=raKey)
        ave_dec = allMatches.aggregate(np.mean, field=decKey)

        # Merge the footprints from the same object together and accumulate
        # information
        object_ids = np.unique(results['object'])
        footprints = []
        all_fluxes = []
        all_flux_errs = []
        num_nobs = []
        self.diaSrcIds = []
        self.diaObjectIds = []

        for id in object_ids:
            mask = results['object'] == id
            num_nobs.append(np.sum(mask))

            footprint = None
            src_ids = []
            fluxes = defaultdict(list)
            flux_errs = defaultdict(list)
            for rec in results[mask]:
                if footprint is None:
                    footprint = rec.getFootprint()
                else:
                    footprint = afwDet.mergeFootprints(footprint, rec.getFootprint())
                src_ids.append(rec.get('id'))

                flux = rec.get(self.config.fluxType)
                if np.isfinite(flux) is False:
                    continue

                filter = self.config.filters[rec.get('filter')]
                calib = self.calibDict[rec.get('visit')][rec.get('ccd')]
                flux_err = rec.get(self.config.fluxType + "Err")
                new_val, new_val_err = scaleFlux(flux, flux_err, calib, self.calib)

                fluxes[filter].append(new_val)
                flux_errs[filter].append(new_val_err)

            self.diaSrcIds.append(src_ids)
            footprints.append(footprint)
            all_fluxes.append(fluxes)
            all_flux_errs.append(flux_errs)

        for i in range(len(ave_ra)):
            rec = cat.addNew()
            self.diaObjectIds.append(rec.get('id'))
            rec.setFootprint(footprints[i])
            rec.set(raKey, ave_ra[i]*geom.radians)
            rec.set(decKey, ave_dec[i]*geom.radians)
            rec.set(nobsKey, num_nobs[i])
            for filter in self.config.filters:
                fluxes = np.array(all_fluxes[i][filter])
                if len(fluxes) == 0:
                    continue

                flux_errs = np.array(all_flux_errs[i][filter])
                flux = self.config.fluxType
                rec.set(keys[f'{flux}_Mean_{filter}'], np.mean(fluxes))
                rec.set(keys[f'{flux}_Sigma_{filter}'], np.std(fluxes, ddof=1))
                rec.set(keys[f'{flux}_Ndata_{filter}'], len(fluxes))
                rec.set(keys[f'{flux}_MeanErr_{filter}'],
                        rec.get(f'{flux}_Sigma_{filter}')/np.sqrt(len(fluxes)))
                residuals = fluxes - rec.get(keys[f'{flux}_Mean_{filter}'])
                rec.set(keys[f'{flux}_Chi2_{filter}'], np.sum((residuals/flux_errs)**2))

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

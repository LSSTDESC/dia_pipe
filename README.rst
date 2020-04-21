########
dia_pipe
########


This package is intended to provide scripts to glue the different pieces of the
LSST difference imaging code into a single pipeline that can be used for DESC
data challenges.  This code has been tested with the LSST stack version w.2018.39.  It
is not compatible with previous releases and has not been tested with later ones.  The
default configuration is designed to work with Data Challenge 1.2 data.

This package requires additional dataset definitions in the obs_lsstCam mapper.  These
are currently availble on the ``dc2/run1.2i_diff`` of obs_lsstCam.

The baseline pipeline includes three steps:
 * Difference image creation at the visit level that produces DIASources and difference images
 * Association of the DIASources into DIAObjects
 * Running forced photometry at the positions of the DIAObjects on all visits


Creating difference images and catalogs
---------------------------------------
It is assumed that you have a repo where the data has been processed all the way through the coadd
and multiBandDriver stages.
Here is the command to create difference images for a set of visits::

    imageDifferenceDriver.py [repo] --rerun [rerun]  --id visit=[visit list]  -C dia_pipe/config/ImageDifferenceDriver.py --cores [cores]

This is a driver script that inherits from ctrl_pool so it can be submitted to a batch system or run
locally with all the standard commands.

For each visit given, it will find the corresponding overlapping coadd tracts/patches.  It is currently
configured to use the science images or ``deepCoadd_calexp`` data products because that is what is available.
This can be changed once there are other coadds produced.  The coadd is then psf-matched to the visit
and subtracted.  The default algorithm is Alard-Lupton, but the Zogy
method is also availible.  More information on the algorithms and settings can be found ??.

The ouptut of this stage is a ``DIASource`` catalog and difference image for every CCD.


Association of DIASources into DIAObjects
-----------------------------------------
The association is done separately for regions on the sky.  For each coadd patch we get all
the overlapping DIASources (within the inner boundary of the tract/patch so there are no duplicate
objects in different catalogs) and associate them into DIAObjects based soley on position.  The association
command is::

    associationDriver.py [repo] --rerun [rerun] --id tract=[tract] filter=u^g^r^i^z^y --cores [cores]

It will fetch DIASources only from the filters explicitly given.  This is also a batch task that will parallelize
over the number of patches in the tract.  Therefore, it doesn't make sense to give more cores than the number
of patches.

There are currently two association algorithms:
 * One based on the MultiMatch algorithm in the LSST stack.  This requires loading all objects into memory before matching and is not 
   meant to handle large number of objects.
 * A simple matching algorithm that keeps a running list of DIAObjects and adds to the list for each new DIASource catalog.  
   This should scale much better than the MultiMatch approach and is currently the default.

It should be fairly simple to add additional association algorithms.

By default the list of visits is chosen by using those that were used to construct the coadd image.
If you want to use a different set of visits you can can append the option::

    --selectId visit=12345^12346^12347


Forced Photometry of DIAObjects
----------------------------------
We can do forced photometry at the positions of the DIAObjects for all the visits::

    forcedPhotCcdDiaDriver.py [repo] --rerun [rerun] --id visit=[visit list] --cores [cores]

This should produce an output for every visit.


Other Forced Photometry Scripts
--------------------------------
There are three other forced photometry algorithms.  The first two will produce an output for each visit.

1. Forced photometry on difference images at the position of the template catalogs.  This will get all template objects 
   from the band the exposure was observed in and within the sky bounding box and perform forced photometry.  Objects that 
   have the position of the parent object outside of the bounding box will be removed from the catalog.  This can be called by::

    forcedPhotCcdTemplateDia [repo] --rerun [rerun] --id visit=[visit list]

   you can also run this via the ``forcedPhotCcdDiaDriver.py`` script by setting a configuration file::

    from lsst.dia.pipe.forcedPhotCcdTemplateDia import ForcedPhotCcdTemplateDiaTask
    config.forced.retarget(ForcedPhotCcdTemplateDiaTask)


2. Forced photometry on the template images at the position of the diaSrc catalog.  This will take 
   the diaSrc catalog and find all template images that contain these objects and perform forced photometry.  
   A ``detect_isPrimary`` flag is added to select unique objects because there will be overlap between different tracts/patches.   This can be called by::

    forcedPhotCoaddDiaSrc [repo] --rerun [rerun] --id visit=[visit list]

   this can also be run via the ``forcedPhotCcdDiaDriver.py`` script by setting a configuration file::

    from lsst.dia.pipe.forcedPhotCoaddDiaSrc import ForcedPhotCoaddDiaSrcTask
    config.forced.retarget(ForcedPhotCoaddDiaSrcTask)


The other forced photometry algorithm performs forced photometry on the template images at the position of the diaObject catalog.  
This is run on tracts and patches and has a different structure than the previous two scripts::

  forcedPhotCoaddDia.py [repo] --rerun [rerun] --id tract=[tract] patch=[patch] filter=[filter]

It has it's own driver script that can be run by::

  forcedPhotCoaddDiaDriver.py [repo] --rerun [rerun] --id tract=[tract] patch=[patch] filter=[filter] --cores [cores]


Setting time and PSF range when running coaddDriver
----------------------------------------------------
This package adds a selector to limit the input visits based on seeing and time that are fed into coaddDriver.py.  Here is an example config file::

  from lsst.dia.pipe.selectImages import TimeBestSeeingWcsSelectImagesTask
  config.select.retarget(TimeBestSeeingWcsSelectImagesTask)
  config.select.minMJD = 59580.0
  config.select.maxMJD = 60310.0
  config.select.minPsfFwhm = 2.5
  config.select.maxPsfFwhm = 3.5

The FWHM values are specified in pixels.



########
dia_pipe
########

.. Add a brief (few sentence) description of what this package provides.

This package is intended to provide scripts to glue the different pieces of the
LSST difference imaging code into a single pipeline that can be used for DESC
data challenges.

The baseline pipeline includes three steps:
* Difference image creation at the visit level that produces DIASources and difference images
* Association of the DIASources into DIAObjects
* Running forced photometry at the positions of the DIAObjects on all visits

It needs to define additional datasets for each camera, so we use a modified version of obs_lsstCam.
The branch is

Creating difference images and catalogs
---------------------------------------
It is assumed that you have a repo where the data has been processed all the way through the coadd 
and multiBandDriver stages.
Here is the command to create difference images for a set of visits:

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
command is:

    associationDriver.py [repo] --rerun [rerun] --id tract=9813 filter=HSC-I^HSC-G^HSC-R^HSC-Z^HSC-Y --cores [cores]

It will fetch DIASources only from the filters explicitly given.  This is also a batch task that will parallelize
over the number of patches in the tract.  Therefore, it doesn't make sense to give more cores than the number
of patches.

There are currently two association algorithms:
* One based on the MultiMatch algorithm in the LSST stack.  This requires loading all objects into memory before
matching and is not meant to handle large number of objects.
* A simple matching algorithm that keeps a running list of DIAObjects and adds to the list for each new DIASource
catalog.  This should scale much better than the MultiMatch approach and is currently the defaul.
It should be fairly simple to add additional association algorithms.


Forced Photometry
----------------------------------
We can now do forced photometry at the positions of the DIAObjects for all the visits.

    forcedPhotCcdDiaDriver.py [repo] --rerun [rerun] --id visit=[visit list] --cores [cores]

This should produce an output for every visit.













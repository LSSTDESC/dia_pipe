Running Difference Imaging at NERSC for Run 1.2p Data
======================================================

This document describes how to run Run1.2p data at NERSC.

Setup software
--------------
The current version of the software is based on the w.2018.39 version of the LSST stack.
I have versions already installed that can be used.::

    source /global/cscratch1/sd/rearmstr/example_diffim/setup.sh
    
    This script will load the appriate modules and setup the lsst software and packages.
    
    module load pe_archive
    module unload python
    module swap PrgEnv-intel PrgEnv-gnu
    module swap gcc gcc/6.3.0
    module rm craype-network-aries
    module rm cray-libsci
    module unload craype
    export CC=gcc
    source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2018_39/loadLSST.bash
    setup lsst_distrib
    setup -jr /global/cscratch1/sd/rearmstr/example_diffim/dia_pipe
    setup -jr /global/cscratch1/sd/rearmstr/example_diffim/obs_lsst
    export OMP_NUM_THREADS=1



Identifying which visits to run
-------------------------------
You need to give the software a list of visits to run.  You can look these up in a database file if you know the tract and patch.
Given that I know one of the DDF tract patches is tract=4849 patch=5,5 I can find the visit list by doing::


    sqlite3 /global/cscratch1/sd/rearmstr/example_diffim/Run1.2_data/rerun/calexp-v4/tracts_mapping.sqlite3

    select distinct(visit) from overlaps where tract=4849 and patch='(5, 5)' order by visit;

This will give you a list of 1074 visits.


Running difference imaging
--------------------------
To enable running on the 1.2p data I had to modify the mapper and registry files to accept a slightly older software version.  This will be updated soon.

To run the difference imaging on one (or more) visits do::

    imageDifferenceDriver.py /global/cscratch1/sd/rearmstr/example_diffim/Run1.2_data/rerun/coadd-v4 --output test_imdiff  --id     visit=431306  -C /global/cscratch1/sd/rearmstr/example_diffim/dia_pipe/config/imageDifferenceDriver.py --batch-type=slurm --mpiexec='-bind-to socket'   --cores 100 --job test --time 500 --batch-options='-C knl -q regular'
    
    the options are:
    - test: is the directory where the output will go
    - cores: specifies the number of cores to use
    - batch-type: if you give "slurm" it will submit to the queue, or "smp" will run on your local machine
    - visit: you can give additional visits by using "^" like "431306^432081..."
    - batch-options: give any other constraints to the batch system
    - time: This is the time per visit.  The total time requested is (number of visits)*time/cores


I have found that splitting up the visits into smaller jobs is much better than trying to submit everything at once.

For every visit, this will get the all the coadd patches that overlap and compute the difference image.


Association
------------
The association of objects is done by looking for all objects that overlap a patch.  You need to give it the tract number and all the filters you want to combine.  It will create a diaObject catalog for each patch.::

    associationDriver.py /global/cscratch1/sd/rearmstr/example_diffim/Run1.2_data/rerun/coadd-v4 --output test_assoc --id tract=4849 filter=u^g^r^i^z^y --cores 10 --batch-type=slurm --mpiexec='-bind-to socket' --cores 100 --job test --time 500 --batch-options='-C knl -q regular'


This program will parallelize over the number of patches/tract.  So it will not be helpful to give it more than the total number of patches which is 81.


Forced photometry
-----------------
This will run forced photometry at the position of the diaObjects.  You must give a list of visits that you want it to run on.::

    forcedPhotCcdDiaDriver.py /global/cscratch1/sd/rearmstr/example_diffim/Run1.2_data/rerun/coadd-v4 --output test --id visit=431306 --cores 10 --clobber-config --clobber-versions --time 100 --batch-options='-C knl -q regular' 




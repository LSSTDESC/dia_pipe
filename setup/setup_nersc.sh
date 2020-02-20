export STACKCVMFS=/cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib
export LSST_STACK_VER=v19.0.0

export DESC_DIA_INSTALL=/global/cfs/cdirs/lsst/groups/SN/dia/code
export DESC_DIA_VER=v19
export DESC_DIA_DIR=$DESC_DIA_INSTALL/$DESC_DIA_VER

module unload python
module swap PrgEnv-intel PrgEnv-gnu
module swap gcc gcc/8.3.0
module rm craype-network-aries
module rm cray-libsci
module unload craype
export CC=gcc
source $STACKCVMFS/$LSST_STACK_VER/loadLSST.bash
setup lsst_distrib
setup -jr $DESC_DIA_DIR/dia_pipe/
setup -jr $DESC_DIA_DIR/obs_lsst/
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1

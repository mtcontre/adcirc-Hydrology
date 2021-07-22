#!/usr/bin/bash

## Load modules
module purge
module load intelc/18.0.0
module load intelfort/18.0.0
module load hdf5/debug_1.10.4_intel-18.0.0
module load mvapich2/2.3b_intel-18.0.0_ch3_ofed-4.1
module load zlib/debug_1.2.11_intel-18.0.0
module list

## Set paths
export PATH=${PATH}:/home/sbunya/opt/netcdf/bin
export NETCDFHOME=`nc-config --prefix`
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NETCDFHOME}/lib
echo 'NETCDFHOME is set to ':
echo "  ${NETCDFHOME}"


#!/usr/bin/bash

## Load modules and set paths
source ./source-me-first.sh

## Make
echo 'Entering cmake'
cmake ../.. -DBUILD_ADCPREP=ON -DBUILD_ADCIRC=ON -DBUILD_ADCSWAN=ON -DBUILD_PADCIRC=ON -DBUILD_PADCSWAN=ON -DENABLE_OUTPUT_NETCDF=ON -DENABLE_MASSCHECK=ON -DDEBUG_ALL_TRACE=ON
echo 'Entering make'
#make -j4 adcprep adcirc adcswan padcirc padcswan
make -j4 adcirc


#!/bin/bash

module2 load gsl/2.3
module2 load openmpi
module2 load gcc/8.1.0
export LD_LIBRARY_PATH=/projects/lanl/Carrington/netcdf/lib:${LD_LIBRARY_PATH}
./Config.pl -install -scheme

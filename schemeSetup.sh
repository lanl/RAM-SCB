#!/bin/bash

module load perl/5.24.1
module load gcc/9.3.0
module load gsl/2.6
module load openmpi/4.0.5
module load netcdf

perl Config.pl -install -ncdf -gsl -compiler=gfortran -openmp -mpi=openmpi
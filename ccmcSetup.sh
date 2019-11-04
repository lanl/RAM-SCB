#!/usr/bin/env bash

###
#  INSTALL SCRIPT FOR LAHAINA01
###

module purge
LD_LIBRARY_PATH=
module load gcc-9.2.0
#TODO: add module load for netcdf here
module load netcdf-fortran-4.5.2/gcc-9.2.0
module load openmpi-3.1.4/gcc-9.2.0
conda activate python37

PKG_CONFIG_PATH=/opt/gnu/gsl-2.5/lib/pkgconfig/
LD_LIBRARY_PATH=/opt/gnu/gsl-2.5/lib:$LD_LIBRARY_PATH
NCDF=/act/netcdf-fortran-4.5.2/gcc-9.2.0

./Config.pl -install -compiler=gfortran -mpi=openmpi -ncdf=${NCDF} -gsl=/opt/gnu/gsl-2.5 -openmp
#now fix the GSL include and link stuff in the make system
sed -i 's/^\(FLAGC.*\)/SEARCH_C = `pkg-config --cflags --libs gsl`\n\1/' Makefile.conf

#Then compile
Make

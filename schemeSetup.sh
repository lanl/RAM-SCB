#!/bin/bash

PSPLINEDIR="${1#*=}"
module2 load netcdf
module2 load netcdf/fortran-4.4.4
module2 load openmpi
echo ${PSPLINEDIR}
./Config.pl -install -compiler=gfortran -mpi=openmpi -pspline=${PSPLINEDIR} -ncdf=/packages2/.packages2/x86_64-pc-linux-gnu-rhel6/netcdf/4.4.4
sed -i '13s/.*/LINK.f90        = ${CUSTOMPATH_MPI}mpif90 -fopenmp/' Makefile.conf
sed -i '34s/.*/OPT3 = -O3 -fopenmp/' Makefile.conf
sed -i '51s/-lngmath/-L\/usr\/lib64\/ncarg\/ -lngmath/' Makefile.conf


#!/bin/bash

module2 load gsl/2.3
module2 load netcdf
module2 load netcdf/fortran-4.4.4
module2 openmpi
./Config.pl -install -scheme
sed -i '66s/.*/.SUFFIXES: .c .f90 .F90 .f .for .ftn .o/' Makefile.conf
sed -i '11iCOMPILE.GCC     = ${CUSTOMPATH_F}gcc' Makefile.conf
sed -i '69i.c.o:' Makefile.conf 
sed -i '70i\\t${COMPILE.f90} ${Cflag3} $<' Makefile.conf
sed -i '14s/.*/LINK.f90        = ${CUSTOMPATH_MPI}mpif90 -fopenmp/' Makefile.conf
sed -i '35s/.*/OPT3 = -O3 -fopenmp/' Makefile.conf

FROM ubuntu:16.04
COPY . /SHIELDS

RUN apt-get update
RUN apt-get install -y git wget python build-essential gfortran libopenmpi-dev openmpi-bin 
RUN apt-get install -y libgsl-dev libgsl2 gsl-bin libgslcblas0 libnetcdf-dev libnetcdff-dev nco netcdf-bin
RUN cd /SHIELDS && ./Config.pl -install -compiler=gfortran -mpi=openmpi -gsl=/usr -ncdf=/usr && sed -i '66s/.*/.SUFFIXES: .c .f90 .F90 .f .for .ftn .o/' Makefile.conf && sed -i '11iCOMPILE.GCC     = ${CUSTOMPATH_F}gcc' Makefile.conf && sed -i '69i.c.o:' Makefile.conf && sed -i '70i\\t${COMPILE.f90} ${Cflag3} $<' Makefile.conf && make

FROM ubuntu:14.04
ENV PERL5LIB=.
COPY . /SHIELDS

RUN apt-get update
RUN apt-get install -y git wget python build-essential gfortran-4.4 gcc-4.4 libopenmpi-dev openmpi-bin 
RUN apt-get install -y libgsl0-dev libgsl0ldbl gsl-bin libgsl0-dbg libnetcdf-dev libnetcdff-dev nco netcdf-bin
RUN cd /SHIELDS && ./Config.pl -install -compiler=gfortran -mpi=openmpi -gsl=/usr -ncdf=/usr && sed -i '66s/.*/.SUFFIXES: .c .f90 .F90 .f .for .ftn .o/' Makefile.conf && sed -i '11iCOMPILE.GCC     = ${CUSTOMPATH_F}gcc' Makefile.conf && sed -i '69i.c.o:' Makefile.conf && sed -i '70i\\t${COMPILE.f90} ${Cflag3} $<' Makefile.conf && make

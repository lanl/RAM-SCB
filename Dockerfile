FROM ubuntu:16.04
ENV PERL5LIB=.
COPY . /SHIELDS

RUN apt-get update
RUN apt-get install -y gfortran 
RUN apt-get install -y git python 
RUN apt-get install -y libopenmpi-dev openmpi-bin
RUN apt-get install -y libgsl-dev libgsl2 gsl-bin libgsl-dbg libnetcdf-dev libnetcdff-dev nco netcdf-bin
RUN gsl-config --version
RUN cd /SHIELDS && ./Config.pl -install -compiler=gfortran -mpi=openmpi -gsl=/usr -ncdf=/usr && sed -i '66s/.*/.SUFFIXES: .c .f90 .F90 .f .for .ftn .o/' Makefile.conf && sed -i '11iCOMPILE.GCC     = ${CUSTOMPATH_F}gcc -std=gnu89' Makefile.conf && sed -i '69i.c.o:' Makefile.conf && sed -i '70i\\t${COMPILE.GCC} ${Cflag3} $<' Makefile.conf && sed -i '35s/.*/OPT3 = -O3 -fopenmp/' Makefile.conf && sed -i '12s/.*/COMPILE.f77 = ${CUSTOMPATH_F}gfortran -std=legacy/' Makefile.conf && sed -i '13s/.*/COMPILE.f90 = ${CUSTOMPATH_F}gfortran -std=legacy/' Makefile.conf && sed -i '14s/.*/LINK.f90 = ${CUSTOMPATH_MPI}mpif90 -fopenmp/' Makefile.conf && make

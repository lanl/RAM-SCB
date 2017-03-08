FROM ubuntu:16.04
COPY . /SHIELDS

RUN apt-get update
RUN apt-get install -y git wget python build-essential gfortran libopenmpi-dev openmpi-bin 
RUN apt-get install -y libnetcdf-dev libnetcdff-dev libncarg-dev
RUN wget http://w3.pppl.gov/rib/repositories/NTCC/files/pspline.tar.gz
RUN mkdir pspline && cd pspline && tar xvf ../pspline.tar.gz &&  make FORTRAN_VARIANT=GCC  NETCDF_DIR=/usr  NETCDF=-lnetcdff
RUN cd /SHIELDS && ./Config.pl -install -compiler=gfortran -mpi=openmpi -pspline=/pspline/LINUX -ncdf=/usr && make



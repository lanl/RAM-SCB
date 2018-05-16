FROM ubuntu:16.04
COPY . /SHIELDS

RUN apt-get update
RUN apt-get install -y git wget python build-essential gfortran libopenmpi-dev openmpi-bin 
RUN apt-get install -y libgsl2 libnetcdf-dev libnetcdff-dev nco netcdf-bin
RUN cd /SHIELDS && ./Config.pl -install -compiler=gfortran -mpi=openmpi -gsl=/usr -ncdf=/usr && make

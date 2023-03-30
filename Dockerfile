FROM ubuntu:20.04
ENV PERL5LIB=.
ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get install -y ssh make gcc git g++ gfortran libopenmpi-dev openmpi-bin libgsl-dev libgsl23 gsl-bin libgsl-dbg libnetcdf-dev libnetcdff-dev nco netcdf-bin

COPY . /SHIELDS
WORKDIR /SHIELDS
RUN ./Config.pl -install -compiler=gfortran -mpi=openmpi -openmp -gsl -ncdf && make

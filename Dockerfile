FROM ubuntu:20.04
ENV PERL5LIB=.
ENV DEBIAN_FRONTEND=noninteractive 
COPY . /SHIELDS

RUN apt-get update
RUN apt-get install -y ssh make git
RUN apt-get install -y gcc g++ gfortran 
RUN apt-get install -y libopenmpi-dev openmpi-bin
RUN apt-get install -y libgsl-dev libgsl23 gsl-bin libgsl-dbg
RUN apt-get install -y libnetcdf-dev libnetcdff-dev nco netcdf-bin

RUN cd /SHIELDS && ./Config.pl -install -compiler=gfortran -mpi=openmpi -openmp -gsl -ncdf && make

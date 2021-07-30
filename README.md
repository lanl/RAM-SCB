# RAM-SCB

[![ramscb-ci](https://github.com/lanl/RAM-SCB/actions/workflows/ramscb-ci.yml/badge.svg)](https://github.com/lanl/RAM-SCB/actions/workflows/ramscb-ci.yml)

The **R**ing current **A**tmosphere interactions **M**odel with **S**elf **C**onsistent magnetic field (**B**) is a unique code that combines a kinetic model of ring current plasma with a three dimensional force-balanced model of the terrestrial magnetic field. The kinetic portion, RAM, solves the kinetic equation to yield the bounce-averaged distribution function as a function of azimuth, radial distance, energy and pitch angle for three ion species (H+, He+, and O+) and, optionally, electrons. The 3-D force balanced magnetic field model, SCB, balances the **J** × **B** force with the divergence of the general pressure tensor to calculate the magnetic field configuration within its domain. The two codes work in tandem, with RAM providing anisotropic pressure to SCB and SCB returning the self-consistent magnetic field through which RAM plasma is advected. RAM-SCB has grown from a research-grade code with limited options and static magnetic field (RAM) to a rich, highly configurable research and operations tool with a multitude of new physics and output products. The RAM-SCB manual provides a guide to users who want to learn how to install, configure, and execute RAM-SCB simulations. While the code is designed to make these steps as straight-forward as possible, it is strongly recommended that users review the publications listed in the Bibliography to ensure a thorough understanding of the physics included in the model.


## Documentation

The RAM-SCB manual with extended installation and usage information can be found at [RAM-SCB/doc/RAM_SCB.pdf](doc/RAM_SCB.pdf).


## Attribution

Researchers who use the RAM-SCB code for scientific research are asked to cite the papers listed below.

1. Jordanova, V. K. et al. (2006), Kinetic simulations of ring current evolution during the Geospace Environment Modeling challenge events, J. Geophys. Res., 111, A11S10, doi:10.1029/2006JA011644.

2. Zaharia, S. et al. (2006), Self-consistent modeling of magnetic fields and plasmas in the inner magnetosphere: Application to the geomagnetic storm, J. Geophys. Res., 111, A11S14, doi:10.1029/2006JA011619.

3. Jordanova, V. K., S. Zaharia, and D. T. Welling (2010), Comparative study of ring current development using empirical, dipolar, and self-consistent magnetic field simulations, J. Geophys. Res., 115(A14):A00J11, doi:10.1029/2010JA015671.

4. Welling, D. T., V. K. Jordanova, S. G. Zaharia, A. Glocer, and G. Toth (2011), The effects of dynamic ionospheric outflow on the ring current, J. Geophys. Res., 116, doi:10.1029/2010JA015642.


## Installation

For full installation details, please see the documentation folder.

To install on Linux (or similar), three main steps are required:
1. Install preqreuisites
  - Fortran/C/C++ compilers are required, we recommend `gfortran`/`gcc`/`g++`
  - GSL is required
  - NetCDF-Fortran is required
  - Perl is required (standard on most linux systems)
2. Configure the build
  - Setup is done using the `Config.pl` Perl script. This sets up the make system and tells RAM-SCB where to find its dependencies.
  - An example setup: `./Config.pl -install -compiler=gfortran -mpi=openmpi -openmp -ncdf -gsl`
  - If `Config.pl` shows the warning "Can't locate share/Scripts/Config.pl in @INC" then add the RAM-SCB directory to the Perl path. In `bash` this is done with ``export PERL5LIB=`pwd` ``
  - Most installations of NetCDF4 and GSL should come with command line utilities to help determine library locations and flags. To auto-detect these use the `-ncdf` and `-gsl` flags without specifying a library location.
3. Compile the code
  - Simply run `make`. If you have multiple cores available for compilation then you can speed things up by running `make -j`


## Usage

```
make rundir RUNDIR=~/desired_run_directory
./ram_scb.exe
```

## Release

This software has been approved for open source release and has been assigned LA-CC-16-077.


## License

The RAM-SCB License can be found at [RAM-SCB/LICENSE.txt](LICENSE.txt).


## Contact

For questions about using RAM-SCB please contact [Vania Jordanova](http://www.lanl.gov/expertise/profiles/view/vania-jordanova).

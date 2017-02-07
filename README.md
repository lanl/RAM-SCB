# RAM-SCB

The **R**ing current **A**tmosphere interactions **M**odel with **S**elf **C**onsistent magnetic field (**B**) is a unique code that combines a kinetic model of ring current plasma with a three dimensional force-balanced model of the terrestrial magnetic field.

## Documentation

The RAM-SCB manual with extended installation and usage information can be found at [RAM-SCB/docs/RAM-SCB.pdf](docs/RAM-SCB.pdf).

## Installation

```
Config.pl -install -compiler=pgf90 -mpi=mpich2
make
````

## Usage

```
make rundir RUNDIR=~/desired_run_directory
./ram_scb.exe
```

## License

The RAM-SCB License can be found at [RAM-SCB/LICENSE.txt](LICENSE.txt).

## Contact

For questions about using RAM-SCB please contact [Vania Jordanova](mailto://vania@lanl.gov).


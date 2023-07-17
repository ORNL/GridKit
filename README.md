# GridKit™

This is experimental code for prototyping interfaces for dynamic 
simulations, sensitivity analysis and optimization. Target applications 
are power grids, but the methodology and the framework could be used 
in other areas without major modifications. 

## Installation Guide

GridKit™ has been built and tested on Linux and Mac platforms. It should
be possible to build it on Windows, as well, with Cygwin or native.
Before installing GridKit™ make sure you have all needed dependencies.

### Dependencies
You should have all of the following installed before installing GridKit™
- A version of
	- [SUNDIALS](https://github.com/LLNL/sundials) >= 6.0.0
	- [Suitesparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) >= 5.x (optional)
	- [Ipopt](https://github.com/coin-or/Ipopt) >= 3.x (optional)
- [CMake](https://cmake.org/) >= 3.12
- C++ 17 compliant compiler

### Installing

GridKit™ uses CMake for build configuration. Per CMake best practices it is recommended 
to build GridKit™ outside the source directory. Building GridKit™ can be as simple as executing
```bash
cmake source_dir
make
make install
```
in the build directory. Dependencies should be autodetected if they are installed in 
standard locations, otherwise you need to specify the location of the dependency 
manually. For example:
```bash
cmake -DSUNDIALS_DIR=/path/to/sundials/install source_dir
```
You can also use `ccmake` or `cmake-gui` tools to adjust GridKit™ build configuration.

### Testing

Several examples are built together with GridKit™ libraries. These are also used
as functionality test and executed by running `ctest` in the build directory.

## Contributors

GridKit™ is written by Slaven Peles (peless@ornl.gov) and has received contributions
from Tamara Becejac (Avangrid), R. Cameron Rutherford (PNNL), Asher J. Mancinelli (NVIDIA), and Reid Gomillion (Virginia Tech).

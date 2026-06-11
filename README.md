# GridKit™

The objective of GridKit™ is to provide a modeling framework for power
systems simulations and analysis that can support multiple advanced
analysys methods, such as dynamic constrained optimization and partitioned
numerical integrators for systems of differential and algebraic equations.
While target applications are power grids, the methodology and the
framework could be used in other areas without major modifications.
GridKit™ supports adding multiple families of models to provide different
representations of power grids and possibly other complex engineered systems.

- Documentation is hosted on [ReadTheDocs](https://gridkit.readthedocs.io/en/latest/)
- Source code is hosted on [GitHub](https://github.com/ORNL/GridKit)

## Installation Guide

GridKit™ has been built and tested on Linux and Mac platforms. It should
be possible to build it on Windows, as well, with Cygwin or native.
Before installing GridKit™ make sure you have all needed dependencies.

### Dependencies
You should have all of the following installed before installing GridKit™
- A version of
    - [SUNDIALS](https://github.com/LLNL/sundials) `develop` branch (optional)
        - To support sparse linear algebra, SUNDIALS must also be built with [KLU support](https://sundials.readthedocs.io/en/latest/sundials/Install_link.html#cmakeoption-ENABLE_KLU). You most likely want this.
    - [Ipopt](https://github.com/coin-or/Ipopt) >= 3.x (optional)
    - [Enzyme](https://github.com/EnzymeAD/Enzyme) >=0.0.206 (optional). Note
      that the use of Enzyme is experimental, and some versions of it have been found to break GridKit code.
        - [LLVM](https://github.com/llvm/llvm-project) >= 15.x. GridKit is
          currently tested with LLVM 16.
    - [ReSolve](https://github.com/ORNL/ReSolve) >= 0.99.3 (optional)
- [CMake](https://cmake.org/) >= 3.13
- C++ 20 compliant compiler

### Cloning and special instructions
```
git clone git@github.com:ORNL/GridKit.git
git submodule update --init third-party/
```
Note, you may need to run the second step periodically as our third party dependencies change.

### Installing

GridKit™ uses CMake for build configuration. The build directory must be
outside the source tree:
```bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ../GridKit
cmake --build .
cmake --install .
```
Dependencies are autodetected if installed in standard locations, otherwise
specify their location explicitly. For example:
```bash
cmake -DSUNDIALS_DIR=/path/to/sundials/install ../GridKit
```
You can also use `ccmake` or `cmake-gui` to adjust the build configuration.

See [INSTALL.md](INSTALL.md) for full instructions, including how to build
each optional dependency and all available CMake options.

### Testing

Several examples are built together with GridKit™ libraries and serve as
functionality tests. Run them with `ctest` in the build directory:
```bash
ctest --output-on-failure
```
To verify the installed CMake configuration works correctly for downstream
projects, run the installation test after `cmake --install`:
```bash
make test_install
```
See [Testing the installation](INSTALL.md#testing-the-installation) in
[INSTALL.md](INSTALL.md) for details.

## Contributors

GridKit™ is written by Slaven Peles (peless@ornl.gov) and has received
contribution from Abdourahman Barry (Virginia Tech), Tamara Becejac (Avangrid),
Adam Birchfield (Texas A&M), Kaleb Brunhoeber (ORNL), Reid Gomillion (Virginia Tech), Nicholson
Koukpaizan (ORNL), Asher J. Mancinelli (NVIDIA), Alex Novotny (Virginia Tech),
Shaked Regev (ORNL), R. Cameron Rutherford (PNNL), and Wiktoria Zielinska (ORNL).

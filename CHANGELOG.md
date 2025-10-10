# Release Changelog

## v0.2

- Added 3, 10, and 37 bus test cases.
- Updated documentation.
- Added JSON parsing.
- Automatic differentiation with enzyme.
- Added PR and issue templates.
- Added Genrou class unit test and example.
- Added Dev container.
- Added back removed tests.
- Added the ability to print matrices to matrix market files for later analysis.
- Added input format specifications.
- Improved IDA interface.
- SUNDIALS interface updates.
- Refactored and reorganized examples.
- Updated variable names and function signatures to follow conventions.
- Added classical generator model.
- Added clang formatting pre-commit.
- CMake fixes.
- Fixed warnings, memory leaks, and failed asserts.
- Added Tgov1 example. 
- Added 10 generator example.
- Improved data structures.
- Removed dead code.
- Added signal node connectivity tools.
- Upgraded to Spack 1.0.0.
- Added `GenClassical` as supported device class in JSON parser.
- Added more context information to JSON parser errors.

## v0.1

- Refactored code to support adding different model families.
- Added PowerElectronics family of models.
- Prototype for sparse system matrix assembly (for PowerElectronics only).
- Added PhasorDynamics family of models.
- Unit testing framework.
- GridKit now depends on Sundials >= 7.0.

## v0.0.7

- Added parser for Matpower files.
- Added system composer, which assembles power flow models from components based on input from Matpower file.
- Updated bus and created generator factory.
- GridKit now depends on Sundials >= 6.0.
- GridKit now requires C++17 compliant compiler.
- Bug fixes.

## v0.0.6

- Refactored CMake to use and export targets
- Bug fixes

## 0.0.5

- Parameter estimation example
- Power flow model composer with test example
- Improved CMake build
- Updated documentation

## 0.0.4

- Add branch component model.
- Add generator with governor component model.
- Add generator connected to infinite bus example.

## 0.0.3

- Add basic system composer.
- Make GridKit forward compatible with Sundials >= 4.0

## 0.0.2

- Add dynamic constraint Ipopt interface and usage example.
- Add Generator 4 component model to the library.
- Require CMake >= 3.0 to build GridKit

## 0.0.1

- Initial release. Skeleton of the library to be developed.

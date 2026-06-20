# Release Changelog

## v0.2

- Added 3, 10, 37, and 39 bus test cases.
- Updated documentation.
- Added support for Read the Docs
- Added JSON parsing.
- Automatic differentiation with enzyme (w.r.t. internal and external variables).
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
- Added `Logger` class.
- Added public function to encapsulate use of nlohmann/json and magic_enum from users.
- Added complete basic consumer example project used for testing installation
- Added `CsrMatrix` class.
- Added `verify` step for signal node links in system model components
- Added support for DependencyTracking::Variable in PowerElectronics models.
- Updated Jacobian value storage from `ScalarT` to `RealT`.
- Added a header file defining constants to be used throughout the code.
- Added `GridKitDocs` target for Doxygen documentation.
- Added a header file defining common math functions (e.g., sigmoid) to be used throughout the code.
- Added capability to print monitored variables in multiple formats, triggered from `Ida::runSimulation`.
- Added IDA statistics object which can be accumulated over multiple simulations.
- Minor performance improvements to residual evaluation in PowerElectronics module.
- Added full support for sparse Jacobians obtained with Enzyme in PhasorDynamics.
- Added `Node` class to the PowerElectronics module to separate nodes from circuit components.
- Refactored Jacobian assembly in `PowerElectronics` module to reuse the CSR pattern.
- Refactored Jacobian assembly in `PhasorDynamics` module to reuse the CSR pattern.
- Removed `COO_Matrix` class use in `PowerElectronics` module.
- Added phasor dynamics application to generalize examples
- Added `LoadZIP` model component type.
- Added component model developer checklist to a README file.
- Added `IEEEST` Stabilizer Model
- Added `SEXS-PTI` Exciter Model
- Added `ESDC1A` Exciter Model
- Added `GENSAL` Machine Model
- Added 200 Bus Synthetic Illinois Case
- Added node objects to `PowerElectronics` module & updated all examples to make use of them.
- Separated internal and external residuals of `PowerElectronics` models.
- Added `CliArgs` class for better management of command-line options.
- Remove data copying between system and components in `PowerElectronics` models.
- Added multi-contingency analysis application.
- Added `BusToSignalAdapter` component for communicating bus voltages and injection currents.
- Added cmake-format hooks, including in pre-commit.
- Added off-nominal tap ratio and phase shift support to the PhasorDynamics `Branch` model.
- Added portable Vector class to GridKit
- Added support for running IDA with fixed time steps

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

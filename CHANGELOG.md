# Release Changelog

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

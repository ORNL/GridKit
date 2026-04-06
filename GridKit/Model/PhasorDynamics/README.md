# Phasor dynamics

This directory contains an implementation of a system model using phasor
dynamics. 

## Modeling checklist

In current GridKit architecture each component model is implemented in six
files:
- `MyModel.cpp`: Basic compilation unit for the case when model does not
                 provide Jacobian.
- `MyModel.hpp`: Component model declaration.
- `MyModelData.hpp`: Definitions of model data including model parameters,
                     model connection ports and model outputs (monitored
                     variables). 
- `MyModelDependencyTracking.cpp`: Compilation unit for Jacobian evaluation
                                   using dependency tracking data type (a.k.a
                                   tapeless automatic differentiation).
- `MyModelEnzyme.cpp`: Compilation unit for Jacobian evaluation with Enzyme.
- `MyModelimpl.hpp`: Actual model implementation.

In essence, a component model developer needs to implement model data
structures in `MyModelData.hpp` file and model equations in `MyModelImpl.hpp`
file. The remaining four files are just supporting infrastructure for
compilation and testing.

We recommend developers follow these steps when adding new component models:
1. Create a subdirectory within appropriate model family directory.
2. Create a README file in markdown format that contains all information
   needed to implement the model. This should include:
    1. List of model parameters in a table format.
    2. List of _derived_ model parameters with mathematical expression
       describing how they are obtained from instantiation parameters.
    3. Model internal variables. Use separate tables for differential and
       algebraic variables.
    4. Model external variables (always algebraic in phasor dynamics).
    5. Model differential and algebraic equations (in separate subsections).
    6. Model initialization procedure with equations in order in which
       initialization computations are performed.
    7. List of model outputs with equations for computing those outputs
       where applicable.
3. Create all six `MyModel*.*pp` implementation files and `CMakeLists.txt`
   file, which specifies build requirements (files to compile, files to
   include, libraries to link and location to install to). Ensure the code
   builds correctly.
    - You may want to start with a "dummy" implementation first to make sure
      the build and installation works correctly before proceeding to the
      implementation.
4. Create unit tests in `tests/UnitTesting/PhasorDynamics` directory. Unit
   test should include:
    1. Constructor tests (smoke tests).
    2. Residual evaluation test (pick input variables to get easy to verify
       residual values).
    3. Jacobian evaluation test (Enzyme generated Jacobian must match Jacobian
       created using dependency tracking variables).
    4. Initialization test (initialize to arbitrary, easy to verify values).
    5. Initialize and evaluate residual test (initialization should ensure
       residual evaluates to zero within prescribed tolerance).
5. Recommended: Create an example in `examples/PhasorDynamics` using the new component.



## Input file parser

Aside from all of the modeling aspects, one additional part of
note is the experimental JSON parsing functionality which can take in a JSON
object (described in `INPUT_FORMAT.md`) and yields a system model. This
is implemented with the [nlohmann/json](https://github.com/nlohmann/json)
parser by adding a `from_json` function for each data structure of note
with the signature `void from_json(const nlohmann::json&, SomeType&)`.

The API reference for the `nlohmann::json` type passed in can be found
[here](https://json.nlohmann.me/api/basic_json) (as well as for the rest
of the library using the navigation on that page). Documentation for some
particular methods of note on this type are here:

- `.at`: https://json.nlohmann.me/api/basic_json/at/
- `.items`: https://json.nlohmann.me/api/basic_json/items/
- `.get_to`: https://json.nlohmann.me/api/basic_json/get_to/

One final thing to note is that the existing parsing code avoids the use
of the `operator[]` implementation on `nlohmann::json` due its lack of
bounds checking. Prefer the use of `.at`.

```cmake
# Optimization solver family: the nonlinear-program contract, dense
# constrained least squares, the Ipopt backend, the rational-approximation
# model space, and the vector fitting estimator. The legacy dynamic
# optimization stubs are retired by this family.
#
# The imported Ipopt target name follows cmake/FindIpopt.cmake.

gridkit_add_library(
  optimization
  SOURCES ConstrainedLeastSquares.cpp IpoptSolver.cpp
  HEADERS OptimizationModel.hpp OptimizationSolver.hpp
          ConstrainedLeastSquares.hpp IpoptSolver.hpp
  LINK_LIBRARIES PRIVATE Ipopt::Ipopt)

add_subdirectory(Rational)
add_subdirectory(VectorFitting)
```

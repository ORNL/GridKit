# Vector fitting: rational approximation posed as constrained optimization.
#
# Not yet wired into GridKit/Solver/CMakeLists.txt; enable together with the
# implementation. Requires GRIDKIT_ENABLE_EIGEN and GRIDKIT_ENABLE_IPOPT.

gridkit_add_library(
  vector_fitting
  SOURCES VectorFitting.cpp
  HEADERS VectorFitting.hpp
  LINK_LIBRARIES PUBLIC GridKit::optimization_rational PRIVATE Eigen3::Eigen)

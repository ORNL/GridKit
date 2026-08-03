```cmake
# Vector fitting: rational approximation posed as constrained optimization.

gridkit_add_library(
  vector_fitting
  SOURCES VectorFitting.cpp
  HEADERS VectorFitting.hpp
  LINK_LIBRARIES PUBLIC GridKit::optimization_rational
                 PRIVATE GridKit::optimization Eigen3::Eigen)
```

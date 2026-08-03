# Rational-approximation model space: sampled responses, rational models,
# delay preprocessing, state-space realization, and passivity assessment.
#
# Not yet wired into GridKit/Solver/CMakeLists.txt; enable together with the
# implementation. Requires GRIDKIT_ENABLE_EIGEN.

gridkit_add_library(
  optimization_rational
  SOURCES RationalModel.cpp MinimumPhase.cpp StateSpaceRealization.cpp
          Passivity.cpp
  HEADERS SampledResponse.hpp RationalModel.hpp MinimumPhase.hpp
          StateSpaceRealization.hpp Passivity.hpp
  LINK_LIBRARIES PRIVATE Eigen3::Eigen)

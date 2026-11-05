gridkit_add_library(
  emt_vector_fit
  SOURCES VectorFitDependencyTracking.cpp
  HEADERS VectorFit.hpp VectorFitData.hpp
  LINK_LIBRARIES PUBLIC GridKit::emt_core)

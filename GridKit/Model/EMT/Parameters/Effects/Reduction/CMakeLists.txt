gridkit_add_library(
  emt_reduction
  SOURCES SeriesReduction.cpp ShuntReduction.cpp
  HEADERS ReductionMap.hpp SeriesReduction.hpp ShuntReduction.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_parameter_views
    GridKit::emt_conductor
    PRIVATE
    Eigen3::Eigen)

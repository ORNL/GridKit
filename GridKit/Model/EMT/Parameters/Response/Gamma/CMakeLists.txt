gridkit_add_library(
  emt_gamma
  SOURCES Gamma.cpp
  HEADERS Gamma.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_parameter_views
    PRIVATE
    Eigen3::Eigen)

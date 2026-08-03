gridkit_add_library(
  emt_zc
  SOURCES Zc.cpp
  HEADERS Zc.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_parameter_views
    PRIVATE
    Eigen3::Eigen)

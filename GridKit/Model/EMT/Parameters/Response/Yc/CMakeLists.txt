gridkit_add_library(
  emt_yc
  SOURCES Yc.cpp
  HEADERS Yc.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_parameter_views
    PRIVATE
    Eigen3::Eigen)

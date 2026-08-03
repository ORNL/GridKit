gridkit_add_library(
  emt_yc
  SOURCES Yc.cpp
  HEADERS Yc.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_series_impedance
    GridKit::emt_shunt_admittance
    PRIVATE
    Eigen3::Eigen)

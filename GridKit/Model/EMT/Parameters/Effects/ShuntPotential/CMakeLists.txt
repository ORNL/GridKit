gridkit_add_library(
  emt_shunt_potential
  SOURCES ShuntPotential.cpp
  HEADERS ShuntPotential.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_tower
    GridKit::emt_conductor)

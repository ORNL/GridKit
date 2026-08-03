gridkit_add_library(
  emt_shunt_admittance
  SOURCES ShuntAdmittance.cpp
  HEADERS ShuntAdmittance.hpp
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_shunt_potential)

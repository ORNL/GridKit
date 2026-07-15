gridkit_add_library(
  power_elec_voltagesource
  SOURCES VoltageSource.cpp
  HEADERS VoltageSource.hpp
  LINK_LIBRARIES GridKit::dense_vector)

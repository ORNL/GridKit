gridkit_add_library(
  power_elec_disgen
  SOURCES DistributedGenerator.cpp
  HEADERS DistributedGenerator.hpp
  LINK_LIBRARIES GridKit::dense_vector)

gridkit_add_library(
  power_elec_tranline
  SOURCES TransmissionLine.cpp
  HEADERS TransmissionLine.hpp
  LINK_LIBRARIES GridKit::dense_vector)

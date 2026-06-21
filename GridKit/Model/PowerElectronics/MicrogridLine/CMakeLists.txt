gridkit_add_library(
  power_elec_microline
  SOURCES MicrogridLine.cpp
  HEADERS MicrogridLine.hpp
  LINK_LIBRARIES GridKit::dense_vector)

gridkit_add_library(
  power_elec_microload
  SOURCES MicrogridLoad.cpp
  HEADERS MicrogridLoad.hpp
  LINK_LIBRARIES GridKit::dense_vector)

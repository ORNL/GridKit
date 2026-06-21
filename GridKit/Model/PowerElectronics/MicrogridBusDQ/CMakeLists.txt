gridkit_add_library(
  power_elec_microbusdq
  SOURCES MicrogridBusDQ.cpp
  HEADERS MicrogridBusDQ.hpp
  LINK_LIBRARIES GridKit::dense_vector)

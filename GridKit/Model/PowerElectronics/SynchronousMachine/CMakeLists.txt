gridkit_add_library(
  power_elec_synmachine
  SOURCES SynchronousMachine.cpp
  HEADERS SynchronousMachine.hpp
  LINK_LIBRARIES GridKit::dense_vector)

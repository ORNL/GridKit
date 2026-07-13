gridkit_add_library(
  power_elec_inductor
  SOURCES Inductor.cpp
  HEADERS Inductor.hpp
  LINK_LIBRARIES GridKit::dense_vector)

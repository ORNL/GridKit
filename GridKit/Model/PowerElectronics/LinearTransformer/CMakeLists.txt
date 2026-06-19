gridkit_add_library(
  power_elec_lineartransformer
  SOURCES LinearTransformer.cpp
  HEADERS LinearTransformer.hpp
  LINK_LIBRARIES GridKit::dense_vector)

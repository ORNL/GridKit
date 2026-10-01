gridkit_add_library(
  power_elec_inductionmotor
  SOURCES InductionMotor.cpp
  HEADERS InductionMotor.hpp
  LINK_LIBRARIES GridKit::dense_vector)

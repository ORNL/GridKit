gridkit_add_library(
  power_elec_bus_line_interface
  SOURCES BusLinePartitionInterface.cpp
  HEADERS BusLinePartitionInterface.hpp PartitionInterface.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::dense_vector
    GridKit::utilities_logger
    GridKit::power_elec_microline)

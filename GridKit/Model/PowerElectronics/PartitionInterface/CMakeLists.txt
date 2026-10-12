gridkit_add_library(
  power_elec_partition_interfaces
  SOURCES BusPartitionInterface.cpp
  HEADERS BusPartitionInterface.hpp PartitionInterface.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::dense_vector
    PUBLIC
    GridKit::utilities_logger)

gridkit_add_library(
  power_elec_bus_line_interface
  SOURCES BusLinePartitionInterface.cpp
  HEADERS BusLinePartitionInterface.hpp
  LINK_LIBRARIES
    PUBLIC
      GridKit::dense_vector
      GridKit::utilities_logger
      GridKit::power_elec_microline
)

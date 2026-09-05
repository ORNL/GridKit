gridkit_add_library(
  power_elec_partition_interfaces
  SOURCES BusPartitionInterface.cpp
  HEADERS BusPartitionInterface.hpp PartitionInterface.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::dense_vector
    PUBLIC
    GridKit::utilities_logger)

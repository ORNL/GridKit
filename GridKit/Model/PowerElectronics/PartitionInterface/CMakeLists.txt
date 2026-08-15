gridkit_add_library(
  power_elec_partition_interfaces
  SOURCES BusPartitionInterface.cpp
  HEADERS BusPartitionInterface.hpp PartitionInterface.hpp
  LINK_LIBRARIES GridKit::dense_vector)

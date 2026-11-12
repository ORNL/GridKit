gridkit_add_library(
  phasor_dynamics_partition
  SOURCES PartitionData.cpp PartitionPlacement.cpp
  HEADERS PartitionData.hpp PartitionPlacement.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_systemmodel
    PUBLIC
    GridKit::utilities_logger
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

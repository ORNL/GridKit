find_package(OpenMP)

add_executable(partition Partition.cpp)

target_link_libraries(
  partition
  PRIVATE
    GridKit::solvers_dyn
    GridKit::power_elec_partition_interfaces
)

add_executable(partitionMicrogrid PartitionMicrogrid.cpp)

target_link_libraries(
  partitionMicrogrid
  PRIVATE
    GridKit::power_elec_disgen
    GridKit::power_elec_microline
    GridKit::power_elec_microload
    GridKit::solvers_dyn
    GridKit::power_elec_microbusdq
    GridKit::power_elec_partition_interfaces
)

add_test(
  NAME PartitionMicrogrid
  COMMAND partitionMicrogrid
)

install(
  TARGETS
    partition
    partitionMicrogrid
  RUNTIME DESTINATION bin
)

if(OpenMP_CXX_FOUND)
  add_executable(
    PartitionScaleMicrogrid
    PartitionScaleMicrogrid.cpp
  )

  target_link_libraries(
    PartitionScaleMicrogrid
    PRIVATE
      GridKit::power_elec_disgen
      GridKit::power_elec_microline
      GridKit::power_elec_microload
      GridKit::solvers_dyn
      GridKit::power_elec_microbusdq
      GridKit::power_elec_partition_interfaces
      OpenMP::OpenMP_CXX
  )

  add_test(
    NAME PartitionScaleMicrogrid
    COMMAND PartitionScaleMicrogrid
  )

  install(
    TARGETS PartitionScaleMicrogrid
    RUNTIME DESTINATION bin
  )
endif()
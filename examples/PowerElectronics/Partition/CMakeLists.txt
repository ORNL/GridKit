add_executable(PartitionMicrogrid PartitionMicrogrid.cpp)
add_executable(PartitionScaleMicrogrid PartitionScaleMicrogrid.cpp)

target_link_libraries(
  PartitionMicrogrid
  PRIVATE GridKit::solvers_dyn GridKit::testing power_elec_partition_examples_helper)

target_link_libraries(
  PartitionScaleMicrogrid
  PRIVATE GridKit::solvers_dyn GridKit::testing power_elec_partition_examples_helper)

if(GRIDKIT_ENABLE_OPENMP)
  target_link_libraries(PartitionMicrogrid PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(PartitionScaleMicrogrid PRIVATE OpenMP::OpenMP_CXX)
endif()

add_test(NAME PartitionMicrogrid COMMAND PartitionMicrogrid)
add_test(NAME PartitionScaleMicrogrid COMMAND PartitionScaleMicrogrid)

install(TARGETS PartitionMicrogrid RUNTIME DESTINATION bin)
install(TARGETS PartitionScaleMicrogrid RUNTIME DESTINATION bin)

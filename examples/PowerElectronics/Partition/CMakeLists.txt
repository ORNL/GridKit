find_package(OpenMP)

add_executable(PartitionMicrogrid PartitionMicrogrid.cpp)
target_include_directories(
  PartitionMicrogrid
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(
  PartitionMicrogrid
  PRIVATE GridKit::power_elec_disgen
          GridKit::power_elec_microline
          GridKit::power_elec_microload
          GridKit::solvers_dyn
          GridKit::power_elec_microbusdq
          GridKit::power_elec_partition_interfaces)

add_executable(PartitionScaleMicrogrid PartitionScaleMicrogrid.cpp)
target_include_directories(
  PartitionScaleMicrogrid
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(
  PartitionScaleMicrogrid
  PRIVATE GridKit::power_elec_disgen
          GridKit::power_elec_microline
          GridKit::power_elec_microload
          GridKit::solvers_dyn
          GridKit::power_elec_microbusdq
          GridKit::power_elec_partition_interfaces)

if(OpenMP_CXX_FOUND)
  target_link_libraries(PartitionScaleMicrogrid PRIVATE OpenMP::OpenMP_CXX)
endif()

add_test(NAME PartitionMicrogrid COMMAND PartitionMicrogrid)
add_test(NAME PartitionScaleMicrogrid COMMAND PartitionScaleMicrogrid)

install(TARGETS PartitionMicrogrid RUNTIME DESTINATION bin)
install(TARGETS PartitionScaleMicrogrid RUNTIME DESTINATION bin)


find_package(OpenMP REQUIRED)

add_executable(partition Partition.cpp)
target_link_libraries(partition     GridKit::solvers_dyn
                                    GridKit::power_elec_partition_interfaces)


add_executable(partitionMicrogrid PartitionMicrogrid.cpp)
target_link_libraries(partitionMicrogrid    GridKit::power_elec_disgen
                                            GridKit::power_elec_microline
                                            GridKit::power_elec_microload
                                            GridKit::solvers_dyn
                                            GridKit::power_elec_microbusdq
                                            GridKit::power_elec_partition_interfaces)


add_executable(PartitionScaleMicrogrid PartitionScaleMicrogrid.cpp)
target_link_libraries(PartitionScaleMicrogrid   GridKit::power_elec_disgen
                                                GridKit::power_elec_microline
                                                GridKit::power_elec_microload
                                                GridKit::solvers_dyn
                                                GridKit::power_elec_microbusdq
                                                GridKit::power_elec_partition_interfaces
                                                OpenMP::OpenMP_CXX)
add_executable(partition Partition.cpp)
target_link_libraries(partition     GridKit::solvers_dyn
                                    GridKit::power_elec_partition_interfaces)


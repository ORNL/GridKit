add_library(
  power_elec_microgrid_network
  INTERFACE
  MicrogridNetwork.hpp
  SystemAssembler.hpp)

target_link_libraries(
  power_elec_microgrid_network
  INTERFACE GridKit::power_elec_disgen
            GridKit::power_elec_microbusdq
            GridKit::power_elec_microline
            GridKit::power_elec_microload)

add_library(
  power_elec_partition_examples_helper
  INTERFACE
  PartitionUtilities.hpp
  JacTestHelper.hpp)

target_link_libraries(
  power_elec_partition_examples_helper
  INTERFACE GridKit::power_elec_partition_interfaces
            GridKit::testing
            GridKit::sparse_matrix
            power_elec_microgrid_network)

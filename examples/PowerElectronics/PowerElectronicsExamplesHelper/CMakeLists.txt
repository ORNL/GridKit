add_library(power_elec_microgrid_network INTERFACE)

add_library(GridKit::power_elec_microgrid_network ALIAS power_elec_microgrid_network)

target_link_libraries(
  power_elec_microgrid_network
  INTERFACE GridKit::power_elec_disgen
            GridKit::power_elec_microbusdq
            GridKit::power_elec_microline
            GridKit::power_elec_microload)
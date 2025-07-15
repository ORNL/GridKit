# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_library(gridkit_phasor_dynamic_components INTERFACE)

add_subdirectory(Branch)
add_subdirectory(Bus)
add_subdirectory(BusFault)
add_subdirectory(Governor)
add_subdirectory(Load)
add_subdirectory(SynchronousMachine)

target_link_libraries(gridkit_phasor_dynamic_components INTERFACE
    GRIDKIT::phasor_dynamics_branch
    GRIDKIT::phasor_dynamics_bus
    GRIDKIT::phasor_dynamics_bus_fault
    GRIDKIT::phasor_dynamics_load
    GRIDKIT::phasor_dynamics_genrou
    GRIDKIT::phasor_dynamics_gen_classical
)
add_library(GRIDKIT::phasor_dynamic_components ALIAS gridkit_phasor_dynamic_components)
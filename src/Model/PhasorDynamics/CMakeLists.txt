# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_library(gridkit_phasor_dynamics_components INTERFACE)
add_library(gridkit_phasor_dynamics_components_dependency_tracking INTERFACE)

add_subdirectory(Branch)
add_subdirectory(Bus)
add_subdirectory(BusFault)
add_subdirectory(Exciter)
add_subdirectory(Governor)
add_subdirectory(Load)
add_subdirectory(SynchronousMachine)

add_subdirectory(SignalNode)

add_library(GRIDKIT::phasor_dynamics_components
  ALIAS gridkit_phasor_dynamics_components)
add_library(GRIDKIT::phasor_dynamics_components_dependency_tracking
  ALIAS gridkit_phasor_dynamics_components_dependency_tracking)

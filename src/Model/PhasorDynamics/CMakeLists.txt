# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_library(phasor_dynamics_components INTERFACE)
add_library(phasor_dynamics_components_dependency_tracking INTERFACE)

add_subdirectory(Branch)
add_subdirectory(Bus)
add_subdirectory(BusFault)
add_subdirectory(Exciter)
add_subdirectory(Governor)
add_subdirectory(Load)
add_subdirectory(SynchronousMachine)

add_subdirectory(SignalNode)

install(
  FILES
    BusBase.hpp
    Component.hpp
    ComponentData.hpp
    ComponentLibrary.hpp
    ComponentSignals.hpp
    SystemModel.hpp
    SystemModelData.hpp
  DESTINATION
    include/GridKit/Model/PhasorDynamics)

add_library(GRIDKIT::phasor_dynamics_components
  ALIAS phasor_dynamics_components)
add_library(GRIDKIT::phasor_dynamics_components_dependency_tracking
  ALIAS phasor_dynamics_components_dependency_tracking)

install(TARGETS
  phasor_dynamics_components
  phasor_dynamics_components_dependency_tracking
  EXPORT gridkit-targets)

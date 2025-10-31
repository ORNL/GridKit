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

gridkit_add_library(phasor_dynamics_core
  SOURCES
    SystemModelData.cpp
  HEADERS
    BusBase.hpp
    Component.hpp
    ComponentData.hpp
    ComponentLibrary.hpp
    ComponentSignals.hpp
    SystemModel.hpp
    SystemModelData.hpp
  INCLUDE_DIRECTORIES
    PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
    PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_core)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_core)

add_library(GridKit::phasor_dynamics_components
  ALIAS phasor_dynamics_components)
add_library(GridKit::phasor_dynamics_components_dependency_tracking
  ALIAS phasor_dynamics_components_dependency_tracking)

install(TARGETS
  phasor_dynamics_components
  phasor_dynamics_components_dependency_tracking
  EXPORT gridkit-targets)

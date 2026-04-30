# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_library(phasor_dynamics_components INTERFACE)
add_library(phasor_dynamics_components_dependency_tracking INTERFACE)

add_library(phasor_dynamics_utils INTERFACE)
target_include_directories(phasor_dynamics_utils INTERFACE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
add_library(GridKit::phasor_dynamics_utils ALIAS phasor_dynamics_utils)

gridkit_add_library(phasor_dynamics_core
  SOURCES
    SystemModelData.cpp
  HEADERS
    BusBase.hpp
    Component.hpp
    ComponentData.hpp
    ComponentLibrary.hpp
    ComponentSignals.hpp
    ConnectedElement.hpp
    GridElement.hpp
    SystemModelData.hpp
  LINK_LIBRARIES
    PUBLIC GridKit::definitions
    PUBLIC GridKit::utilities_logger
    PUBLIC GridKit::sparse_matrix
    PRIVATE GridKit::phasor_dynamics_utils)

target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_core)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_core)

gridkit_add_library(phasor_dynamics_system_model
  SOURCES
    SystemModel.cpp
  HEADERS
    SystemModel.hpp
  LINK_LIBRARIES
    PUBLIC GridKit::phasor_dynamics_core
    PRIVATE GridKit::phasor_dynamics_utils)

gridkit_add_library(phasor_dynamics_system_model_dependency_tracking
  SOURCES
    SystemModelDependencyTracking.cpp
  LINK_LIBRARIES
    PUBLIC GridKit::phasor_dynamics_core
  PRIVATE GridKit::phasor_dynamics_utils)

target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_system_model)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_system_model_dependency_tracking)

add_subdirectory(Branch)
add_subdirectory(Bus)
add_subdirectory(BusFault)
add_subdirectory(Exciter)
add_subdirectory(Governor)
add_subdirectory(Load)
add_subdirectory(LoadZIP)
add_subdirectory(Stabilizer)
add_subdirectory(SynchronousMachine)

add_subdirectory(SignalNode)

add_library(GridKit::phasor_dynamics_components
  ALIAS phasor_dynamics_components)
add_library(GridKit::phasor_dynamics_components_dependency_tracking
  ALIAS phasor_dynamics_components_dependency_tracking)

install(TARGETS
  phasor_dynamics_components
  phasor_dynamics_components_dependency_tracking
  EXPORT gridkit-targets)

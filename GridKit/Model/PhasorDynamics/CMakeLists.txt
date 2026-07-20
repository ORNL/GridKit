# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_library(phasor_dynamics_components INTERFACE)
add_library(phasor_dynamics_components_dependency_tracking INTERFACE)

gridkit_add_library(
  phasor_dynamics_core
  INTERFACE_TARGET
  HEADERS BusBase.hpp
          Component.hpp
          ComponentData.hpp
          ComponentLibrary.hpp
          ComponentSignals.hpp
  LINK_LIBRARIES
    INTERFACE
    GridKit::definitions
    INTERFACE
    GridKit::utilities_logger
    INTERFACE
    GridKit::dense_vector
    INTERFACE
    GridKit::sparse_matrix)

target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_core)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_core)

add_subdirectory(Branch)
add_subdirectory(Bus)
add_subdirectory(BusFault)
add_subdirectory(BusToSignalAdapter)
add_subdirectory(Exciter)
add_subdirectory(Governor)
add_subdirectory(Load)
add_subdirectory(SignalNode)
add_subdirectory(SignalSource)
add_subdirectory(Stabilizer)
add_subdirectory(SynchronousMachine)

add_library(GridKit::phasor_dynamics_components ALIAS phasor_dynamics_components)
add_library(GridKit::phasor_dynamics_components_dependency_tracking ALIAS
            phasor_dynamics_components_dependency_tracking)


if(GRIDKIT_ENABLE_ENZYME)
  set(_system_model_source SystemModelEnzyme.cpp)
else()
  set(_system_model_source SystemModel.cpp)
endif()

gridkit_add_library(
  phasor_dynamics_systemmodel
  SOURCES ${_system_model_source} SystemModelData.cpp
  HEADERS SystemModel.hpp SystemModelData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_components
  INCLUDE_DIRECTORIES
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

gridkit_add_library(
  phasor_dynamics_systemmodel_dependency_tracking
  SOURCES SystemModelDependencyTracking.cpp
  HEADERS SystemModel.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_components_dependency_tracking)

install(TARGETS phasor_dynamics_components phasor_dynamics_components_dependency_tracking
        EXPORT gridkit-targets)

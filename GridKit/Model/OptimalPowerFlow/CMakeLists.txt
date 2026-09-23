add_library(optimal_power_flow_components INTERFACE)
add_library(optimal_power_flow_components_dependency_tracking INTERFACE)

gridkit_add_library(
  optimal_power_flow_core
  INTERFACE_TARGET
  HEADERS Component.hpp ComponentData.hpp ComponentModel.hpp
  LINK_LIBRARIES
    INTERFACE
    GridKit::definitions
    INTERFACE
    GridKit::Utilities
    INTERFACE
    GridKit::utilities_logger
    INTERFACE
    GridKit::model_state_data)

target_link_libraries(
  optimal_power_flow_components
  INTERFACE GridKit::optimal_power_flow_core)
target_link_libraries(
  optimal_power_flow_components_dependency_tracking
  INTERFACE GridKit::optimal_power_flow_core)

add_subdirectory(Branch)
add_subdirectory(Bus)
add_subdirectory(Generator)
add_subdirectory(Load)
add_subdirectory(Shunt)

add_library(GridKit::optimal_power_flow_components ALIAS optimal_power_flow_components)
add_library(GridKit::optimal_power_flow_components_dependency_tracking ALIAS
            optimal_power_flow_components_dependency_tracking)

gridkit_add_library(
  optimal_power_flow_systemmodel
  SOURCES SystemModelEnzyme.cpp SystemModelData.cpp MatpowerData.cpp
  HEADERS SystemModel.hpp SystemModelData.hpp MatpowerData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::optimal_power_flow_components
    PUBLIC
    GridKit::dense_vector
    PUBLIC
    GridKit::sparse_matrix)

gridkit_add_library(
  optimal_power_flow_systemmodel_dependency_tracking
  SOURCES SystemModelDependencyTracking.cpp
  HEADERS SystemModel.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::optimal_power_flow_components_dependency_tracking
    PUBLIC
    GridKit::dense_vector
    PUBLIC
    GridKit::sparse_matrix)

install(TARGETS optimal_power_flow_components optimal_power_flow_components_dependency_tracking
        EXPORT gridkit-targets)

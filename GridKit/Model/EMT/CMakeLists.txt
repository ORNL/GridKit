# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#]]

add_library(emt_components INTERFACE)
add_library(emt_components_dependency_tracking INTERFACE)

gridkit_add_library(
  emt_core
  INTERFACE_TARGET
  HEADERS Component.hpp
          ComponentData.hpp
          ComponentLibrary.hpp
          ComponentSignals.hpp
          PhasorInitialization.hpp
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
  emt_components
  INTERFACE GridKit::emt_core)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_core)

add_subdirectory(Signal)
add_subdirectory(Operators)
add_subdirectory(Component)

add_library(GridKit::emt_components ALIAS emt_components)
add_library(GridKit::emt_components_dependency_tracking ALIAS emt_components_dependency_tracking)

if(GRIDKIT_ENABLE_ENZYME)
  set(_system_model_source SystemModelEnzyme.cpp)
else()
  set(_system_model_source SystemModel.cpp)
endif()

gridkit_add_library(
  emt_systemmodel
  SOURCES ${_system_model_source} SystemModelData.cpp StateDataAdapter.cpp
  HEADERS Container.hpp
          ContainerData.hpp
          SystemModel.hpp
          SystemModelData.hpp
          StateDataAdapter.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_components
    PUBLIC
    GridKit::model_state_data
  INCLUDE_DIRECTORIES
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

gridkit_add_library(
  emt_systemmodel_dependency_tracking
  SOURCES SystemModelDependencyTracking.cpp
  HEADERS SystemModel.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_components_dependency_tracking
  INCLUDE_DIRECTORIES
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

install(TARGETS emt_components emt_components_dependency_tracking EXPORT gridkit-targets)

# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

add_library(emt_components INTERFACE)
add_library(emt_components_dependency_tracking INTERFACE)

gridkit_add_library(
  emt_dae_analysis
  SOURCES DaeAnalysis.cpp
  HEADERS DaeAnalysis.hpp)
if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  target_link_libraries(emt_dae_analysis PRIVATE SUNDIALS::KLU)
  target_compile_definitions(emt_dae_analysis PRIVATE GRIDKIT_ENABLE_SUNDIALS_SPARSE)
endif()

gridkit_add_library(
  emt_core
  INTERFACE_TARGET
  HEADERS Component.hpp
          ComponentData.hpp
          ComponentInitialization.hpp
          ComponentLibrary.hpp
          ComponentSignals.hpp
          InitialState.hpp
          JsonValidation.hpp
          PhasorInitialization.hpp
  INCLUDE_DIRECTORIES INTERFACE $<BUILD_INTERFACE:${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include>
                      $<BUILD_INTERFACE:${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include>
  LINK_LIBRARIES
    INTERFACE
    GridKit::definitions
    INTERFACE
    GridKit::utilities_logger
    INTERFACE
    GridKit::dense_vector
    INTERFACE
    GridKit::sparse_matrix)

install(DIRECTORY ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include/magic_enum DESTINATION include)
install(FILES ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/LICENSE
        DESTINATION share/licenses/GridKit/magic-enum)
install(DIRECTORY ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include/nlohmann DESTINATION include)
install(FILES ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/LICENSE.MIT
        DESTINATION share/licenses/GridKit/nlohmann-json)

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
  SOURCES ${_system_model_source} SystemModelData.cpp
  HEADERS Container.hpp
          ContainerRuntime.hpp
          ContainerData.hpp
          SystemModel.hpp
          SystemModelData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_components
    PUBLIC
    GridKit::emt_dae_analysis
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
    PUBLIC
    GridKit::emt_dae_analysis
  INCLUDE_DIRECTORIES
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

install(TARGETS emt_components emt_components_dependency_tracking EXPORT gridkit-targets)

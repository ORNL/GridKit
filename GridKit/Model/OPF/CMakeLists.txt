gridkit_add_library(
  model_opf_core
  INTERFACE_TARGET
  HEADERS Component.hpp
  LINK_LIBRARIES
    INTERFACE
    GridKit::model_state_data
    INTERFACE
    GridKit::optimization_core)

gridkit_add_library(model_opf_components INTERFACE_TARGET)

add_subdirectory(Bus)
add_subdirectory(Generator)
add_subdirectory(Load)
add_subdirectory(Shunt)
add_subdirectory(Branch)

gridkit_add_library(
  model_opf
  SOURCES System.cpp SystemData.cpp
  HEADERS System.hpp SystemData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::model_opf_components
    PUBLIC
    GridKit::model_opf_core
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

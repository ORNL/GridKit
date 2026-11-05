if(GRIDKIT_ENABLE_ENZYME)
  add_library(emt_components INTERFACE)

  gridkit_add_library(
    emt_core
    INTERFACE_TARGET
    HEADERS Component.hpp ComponentData.hpp
    LINK_LIBRARIES
      INTERFACE
      GridKit::definitions
      INTERFACE
      GridKit::Utilities
      INTERFACE
      GridKit::dense_vector
      INTERFACE
      GridKit::sparse_matrix
      INTERFACE
      GridKit::utilities_logger)

  target_link_libraries(emt_components INTERFACE GridKit::emt_core)

  add_subdirectory(Bus)
  add_subdirectory(Component)
  add_subdirectory(Operators)

  add_library(GridKit::emt_components ALIAS emt_components)

  gridkit_add_library(
    emt_systemmodel
    SOURCES SystemModel.cpp SystemModelData.cpp
    HEADERS SystemModel.hpp
            SystemModelData.hpp
            ComponentLibrary.hpp
            Monitor.hpp
            StateJSONParser.hpp
    LINK_LIBRARIES
      PUBLIC
      GridKit::emt_components
      PUBLIC
      GridKit::dense_vector
      PUBLIC
      GridKit::sparse_matrix
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

  install(TARGETS emt_components EXPORT gridkit-targets)
endif()

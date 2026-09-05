gridkit_add_library(
  emt_converter
  SOURCES Converter.cpp
  HEADERS Converter.hpp ConverterData.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal)

gridkit_add_library(
  emt_converter_dependency_tracking
  SOURCES ConverterDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_converter)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_converter_dependency_tracking)

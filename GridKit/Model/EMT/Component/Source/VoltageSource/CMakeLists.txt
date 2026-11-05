gridkit_add_library(
  emt_voltage_source
  SOURCES VoltageSourceDependencyTracking.cpp
  HEADERS VoltageSource.hpp VoltageSourceData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_vector_fit)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_voltage_source)

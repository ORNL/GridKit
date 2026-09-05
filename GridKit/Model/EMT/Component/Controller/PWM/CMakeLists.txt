gridkit_add_library(
  emt_pwm
  SOURCES Pwm.cpp
  HEADERS Pwm.hpp PwmData.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal)

gridkit_add_library(
  emt_pwm_dependency_tracking
  SOURCES PwmDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_pwm)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_pwm_dependency_tracking)

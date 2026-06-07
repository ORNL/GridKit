set(_install_headers SampledSignalSource.hpp SampledSignalSourceData.hpp)

gridkit_add_library(
  phasor_dynamics_sampled_signal_source
  SOURCES SampledSignalSource.cpp
  HEADERS ${_install_headers}
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)

gridkit_add_library(
  phasor_dynamics_sampled_signal_source_dependency_tracking
  SOURCES SampledSignalSourceDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_signal_dependency_tracking)

target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_sampled_signal_source)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_sampled_signal_source_dependency_tracking)

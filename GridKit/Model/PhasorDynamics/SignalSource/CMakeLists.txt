set(_install_headers
    ConstantSignalSource.hpp
    ConstantSignalSourceData.hpp
    ForcedOscillation.hpp
    ForcedOscillationData.hpp)

gridkit_add_library(
  phasor_dynamics_signalsource
  SOURCES ConstantSignalSource.cpp ForcedOscillation.cpp
  HEADERS ${_install_headers}
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)

gridkit_add_library(
  phasor_dynamics_signalsource_dependency_tracking
  SOURCES ConstantSignalSourceDependencyTracking.cpp ForcedOscillationDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_signalsource)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_signalsource_dependency_tracking)

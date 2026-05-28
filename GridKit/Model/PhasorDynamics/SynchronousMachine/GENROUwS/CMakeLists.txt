set(_install_headers Genrou.hpp GenrouData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    phasor_dynamics_genrou
    SOURCES GenrouEnzyme.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::phasor_dynamics_core
      PUBLIC
      GridKit::phasor_dynamics_signal
      PRIVATE
      ClangEnzymeFlags
    COMPILE_OPTIONS
      PRIVATE
      -mllvm
      -enzyme-auto-sparsity=1
      -fno-math-errno
      -fno-vectorize)
else()
  gridkit_add_library(
    phasor_dynamics_genrou
    SOURCES Genrou.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::phasor_dynamics_core
      PUBLIC
      GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(
  phasor_dynamics_genrou_dependency_tracking
  SOURCES GenrouDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_genrou)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_genrou_dependency_tracking)

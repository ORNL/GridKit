set(_install_headers
  Gensal.hpp
  GensalData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_gensal
    SOURCES
      GensalEnzyme.cpp
    HEADERS
      ${_install_headers}
    INCLUDE_DIRECTORIES
      PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC GridKit::phasor_dynamics_core
      PUBLIC GridKit::phasor_dynamics_signal
      PRIVATE ClangEnzymeFlags
    COMPILE_OPTIONS
      PRIVATE -mllvm -enzyme-auto-sparsity=1 -fno-math-errno -fno-vectorize)
else()
  gridkit_add_library(phasor_dynamics_gensal
    SOURCES
      Gensal.cpp
    HEADERS
      ${_install_headers}
    INCLUDE_DIRECTORIES
      PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC GridKit::phasor_dynamics_core
      PUBLIC GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(phasor_dynamics_gensal_dependency_tracking
  SOURCES
    GensalDependencyTracking.cpp
  INCLUDE_DIRECTORIES
    PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC GridKit::phasor_dynamics_core
    PUBLIC GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_gensal)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_gensal_dependency_tracking)

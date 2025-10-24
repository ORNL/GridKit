
set(_install_headers
  Genrou.hpp
  GenrouData.hpp)

# Disabled GenrouEnzyme because of an intermittent build issue
#if(GRIDKIT_ENABLE_ENZYME)
#  gridkit_add_library(phasor_dynamics_genrou
#    SOURCES
#      GenrouEnzyme.cpp
#    HEADERS
#      ${_install_headers}
#    LINK_LIBRARIES
#      PUBLIC GRIDKIT::phasor_dynamics_signal ClangEnzymeFlags
#    COMPILE_OPTIONS
#      -mllvm -enzyme-auto-sparsity=1 -fno-math-errno -fno-vectorize)
#else()
  gridkit_add_library(phasor_dynamics_genrou
    SOURCES
      Genrou.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PUBLIC GRIDKIT::phasor_dynamics_signal)
#endif()

gridkit_add_library(phasor_dynamics_genrou_dependency_tracking
  SOURCES
    GenrouDependencyTracking.cpp
  LINK_LIBRARIES
    PUBLIC GRIDKIT::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components
  INTERFACE GRIDKIT::phasor_dynamics_genrou)
target_link_libraries(gridkit_phasor_dynamics_components_dependency_tracking
  INTERFACE GRIDKIT::phasor_dynamics_genrou_dependency_tracking)

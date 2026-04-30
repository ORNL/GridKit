# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

set(_install_headers
  GenClassical.hpp
  GenClassicalData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_gen_classical
    SOURCES
      GenClassicalEnzyme.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PRIVATE ClangEnzymeFlags
      PUBLIC GridKit::phasor_dynamics_core
      PRIVATE GridKit::phasor_dynamics_utils
    COMPILE_OPTIONS
      PRIVATE -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_gen_classical
    SOURCES
      GenClassical.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PRIVATE GridKit::phasor_dynamics_utils
      PUBLIC GridKit::phasor_dynamics_core)
endif()

gridkit_add_library(phasor_dynamics_gen_classical_dependency_tracking
  SOURCES 
    GenClassicalDependencyTracking.cpp
  LINK_LIBRARIES
    PRIVATE GridKit::phasor_dynamics_utils
    PUBLIC GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_gen_classical)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_gen_classical_dependency_tracking)

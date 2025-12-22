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
    INCLUDE_DIRECTORIES
      PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PRIVATE ClangEnzymeFlags
    COMPILE_OPTIONS
      PRIVATE -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_gen_classical
    SOURCES
      GenClassical.cpp
    HEADERS
      ${_install_headers}
    INCLUDE_DIRECTORIES
      PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
endif()

gridkit_add_library(phasor_dynamics_gen_classical_dependency_tracking
  SOURCES GenClassicalDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_gen_classical)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_gen_classical_dependency_tracking)

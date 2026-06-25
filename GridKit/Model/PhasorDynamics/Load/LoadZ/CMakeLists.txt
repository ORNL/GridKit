# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
# ]]

set(_install_headers LoadZ.hpp LoadZData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    phasor_dynamics_loadz
    SOURCES LoadZEnzyme.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PRIVATE
      ClangEnzymeFlags
      PUBLIC
      GridKit::phasor_dynamics_core
    COMPILE_OPTIONS
      PRIVATE
      -mllvm
      -enzyme-auto-sparsity=1
      -fno-math-errno)
else()
  gridkit_add_library(
    phasor_dynamics_loadz
    SOURCES LoadZ.cpp
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core
    HEADERS ${_install_headers})
endif()

gridkit_add_library(
  phasor_dynamics_loadz_dependency_tracking
  SOURCES LoadZDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_loadz)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_loadz_dependency_tracking)



# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
# ]]

set(_install_headers
  Load.hpp
  LoadData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_load
    SOURCES
      LoadEnzyme.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PRIVATE ClangEnzymeFlags
      PUBLIC GridKit::phasor_dynamics_core
      PRIVATE GridKit::phasor_dynamics_utils
    COMPILE_OPTIONS
      PRIVATE -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_load
    SOURCES
      Load.cpp
    LINK_LIBRARIES
      PUBLIC GridKit::phasor_dynamics_core
      PRIVATE GridKit::phasor_dynamics_utils
    HEADERS
      ${_install_headers})
endif()

gridkit_add_library(phasor_dynamics_load_dependency_tracking
  SOURCES 
    LoadDependencyTracking.cpp
  LINK_LIBRARIES
    PRIVATE GridKit::phasor_dynamics_utils
    PUBLIC GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_load)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_load_dependency_tracking)

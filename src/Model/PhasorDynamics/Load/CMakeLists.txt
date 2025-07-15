

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
# ]]


if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_load
    SOURCES
      LoadEnzyme.cpp
    LINK_LIBRARIES
      ClangEnzymeFlags
    COMPILE_OPTIONS
      -mllvm -enzyme-auto-sparsity=1 -fno-math-errno
    OUTPUT_NAME
      gridkit_phasor_dynamics_load)
else()
  gridkit_add_library(phasor_dynamics_load
    SOURCES
      Load.cpp
    OUTPUT_NAME
      gridkit_phasor_dynamics_load)
endif()

gridkit_add_library(phasor_dynamics_load_dependency_tracking
  SOURCES
    LoadDependencyTracking.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_load_dependency_tracking)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components INTERFACE GRIDKIT::phasor_dynamics_load)

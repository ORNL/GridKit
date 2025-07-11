# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_gen_classical
    SOURCES
      GenClassicalEnzyme.cpp
    LINK_LIBRARIES
      ClangEnzymeFlags
    COMPILE_OPTIONS
      -mllvm -enzyme-auto-sparsity=1 -fno-math-errno
    OUTPUT_NAME
      gridkit_phasor_dynamics_gen_classical)
else()
  gridkit_add_library(phasor_dynamics_gen_classical
    SOURCES
      GenClassical.cpp
    OUTPUT_NAME
      gridkit_phasor_dynamics_gen_classical)
endif()

  gridkit_add_library(phasor_dynamics_gen_classical_dependency_tracking
    SOURCES
      GenClassicalDependencyTracking.cpp
    OUTPUT_NAME
      gridkit_phasor_dynamics_gen_classical_dependency_tracking)

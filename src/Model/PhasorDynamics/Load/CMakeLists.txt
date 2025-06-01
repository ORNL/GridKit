 

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
# ]]

gridkit_add_library(phasor_dynamics_load
  SOURCES
    Load.cpp
  LINK_LIBRARIES
    ClangEnzymeFlags
  OUTPUT_NAME
    gridkit_phasor_dynamics_load)


# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

gridkit_add_library(phasor_dynamics_synchronous_machine
  SOURCES
    SynchronousMachine.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_synchronous_machine)

add_subdirectory(GENROUwS)

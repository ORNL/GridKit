# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#    - Adam Birchfield <abirchfield@tamu.edu>
# ]]

gridkit_add_library(phasor_dynamics_governortgov1
  SOURCES
    Tgov1.cpp
  LINK_LIBRARIES
    GRIDKIT::phasor_dynamics_genrou
  OUTPUT_NAME
    gridkit_phasor_dynamics_governortgov1)

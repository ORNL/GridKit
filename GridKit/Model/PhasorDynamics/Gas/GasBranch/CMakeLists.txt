# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers GasBranch.hpp GasBranchData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  set(_gas_branch_sources GasBranchEnzyme.cpp)
else()
  set(_gas_branch_sources GasBranch.cpp)
endif()

gridkit_add_library(
  phasor_dynamics_gas_branch
  SOURCES ${_gas_branch_sources}
  HEADERS ${_install_headers}
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

gridkit_add_library(
  phasor_dynamics_gas_branch_dependency_tracking
  SOURCES GasBranchDependencyTracking.cpp
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_gas_branch)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_gas_branch_dependency_tracking)

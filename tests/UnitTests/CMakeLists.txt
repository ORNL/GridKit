#
add_subdirectory(AutomaticDifferentiation)
add_subdirectory(LinearAlgebra)
add_subdirectory(Math)
add_subdirectory(Network)
add_subdirectory(PhasorDynamics)
add_subdirectory(PowerElectronics)
add_subdirectory(Solver)
add_subdirectory(Utilities)

# Gated behind GridKit_ENABLE_NOTEBOOKS: the pytest suite imports the
# application/notebooks packages, which requires their (non-C++) Python
# dependencies to be installed — see application/notebooks/CMakeLists.txt.
if(GridKit_ENABLE_NOTEBOOKS)
  add_subdirectory(Python)
endif()

#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

set(solvers_opt_sources IpoptSolver.cpp)
set(solvers_opt_headers IpoptSolver.hpp OptimizationResult.hpp OptimizationSolver.hpp)
set(solvers_opt_libraries
    PUBLIC
    IPOPT
    PUBLIC
    GridKit::dense_vector
    PUBLIC
    GridKit::sparse_matrix)

# The evaluator-based Ipopt solver is independent of SUNDIALS. Keep the
# existing dynamic optimization interfaces available when IDAS is enabled.
if(TARGET solvers_dyn)
  list(
    APPEND
    solvers_opt_sources
    DynamicObjective.cpp
    DynamicConstraint.cpp)
  list(
    APPEND
    solvers_opt_headers
    DynamicConstraint.hpp
    DynamicObjective.hpp)
  list(
    APPEND
    solvers_opt_libraries
    PUBLIC
    GridKit::solvers_dyn)
endif()

gridkit_add_library(
  solvers_opt
  SOURCES ${solvers_opt_sources}
  HEADERS ${solvers_opt_headers}
  LINK_LIBRARIES ${solvers_opt_libraries})

# An executable's RUNPATH is not used to locate indirect dependencies. Keep
# the discovered Ipopt location on the installed solver library itself.
if(IPOPT_LIBRARY_DIR)
  set_property(
    TARGET solvers_opt
    APPEND
    PROPERTY INSTALL_RPATH "${IPOPT_LIBRARY_DIR}")
endif()

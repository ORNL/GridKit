#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  solvers_steady
  SOURCES Kinsol.cpp
  HEADERS Kinsol.hpp SteadyStateSolver.hpp
  LINK_LIBRARIES PUBLIC SUNDIALS::kinsol GridKit::dense_vector)

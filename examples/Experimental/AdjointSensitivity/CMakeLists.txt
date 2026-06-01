#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_executable(adjoint AdjointSensitivity.cpp)
target_link_libraries(
  adjoint
  GridKit::bus
  GridKit::generator4
  GridKit::solvers_dyn)
install(TARGETS adjoint DESTINATION bin)

add_test(NAME AdjointSens COMMAND $<TARGET_FILE:adjoint>)

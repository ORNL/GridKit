#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_executable(adjoint AdjointSensitivity.cpp)
target_link_libraries(adjoint GRIDKIT::bus GRIDKIT::generator4 GRIDKIT::solvers_dyn)
install(TARGETS adjoint DESTINATION bin)

add_test(NAME AdjointSens COMMAND $<TARGET_FILE:adjoint> )

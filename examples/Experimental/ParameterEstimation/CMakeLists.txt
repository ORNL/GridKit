#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(paramest ParameterEstimation.cpp)
target_link_libraries(paramest
  GridKit::bus
  GridKit::generator4param
  GridKit::solvers_opt
  GridKit::solvers_dyn
  GridKit::Utilities
)
# Fortran linker needed to link to HSL solvers
set_property(TARGET paramest PROPERTY LINKER_LANGUAGE Fortran)
install(TARGETS paramest RUNTIME DESTINATION bin)

add_test(NAME ParameterEst COMMAND $<TARGET_FILE:paramest>)

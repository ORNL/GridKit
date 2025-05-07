#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(paramest ParameterEstimation.cpp)
target_link_libraries(paramest
  GRIDKIT::bus
  GRIDKIT::generator4param
  GRIDKIT::solvers_opt
  GRIDKIT::solvers_dyn
  GRIDKIT::Utilities
)
# Fortran linker needed to link to HSL solvers
set_property(TARGET paramest PROPERTY LINKER_LANGUAGE Fortran)
install(TARGETS paramest RUNTIME DESTINATION bin)

add_test(NAME ParameterEst COMMAND $<TARGET_FILE:paramest>)

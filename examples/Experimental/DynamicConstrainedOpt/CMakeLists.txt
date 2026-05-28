#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(dynconopt DynamicConstrainedOpt.cpp)
target_link_libraries(
  dynconopt
  GridKit::generator4
  GridKit::generator2
  GridKit::bus
  GridKit::solvers_dyn
  GridKit::solvers_opt)
# Fortran linker needed to link to HSL solvers
set_property(TARGET dynconopt PROPERTY LINKER_LANGUAGE Fortran)
install(TARGETS dynconopt DESTINATION bin)

add_test(NAME DynamicConOpt COMMAND $<TARGET_FILE:dynconopt>)

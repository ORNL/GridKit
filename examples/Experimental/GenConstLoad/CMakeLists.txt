#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(genconstload GenConstLoad.cpp)
target_link_libraries(
  genconstload
  GridKit::generator4governor
  GridKit::bus
  GridKit::load
  GridKit::solvers_dyn
  GridKit::solvers_opt)
# Fortran linker needed to link to HSL solvers
set_property(TARGET genconstload PROPERTY LINKER_LANGUAGE Fortran)
install(TARGETS genconstload DESTINATION bin)

add_test(NAME GenConstLoad COMMAND $<TARGET_FILE:genconstload>)

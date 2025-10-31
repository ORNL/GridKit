#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(geninfbus GenInfiniteBus.cpp)
target_link_libraries(geninfbus
    GridKit::bus
    GridKit::generator4
    GridKit::solvers_opt
    GridKit::solvers_dyn
)
# Fortran linker needed to link to HSL solvers
set_property(TARGET geninfbus PROPERTY LINKER_LANGUAGE Fortran)
install(TARGETS geninfbus DESTINATION bin)

add_test(NAME GenInfiniteBus COMMAND $<TARGET_FILE:geninfbus>)

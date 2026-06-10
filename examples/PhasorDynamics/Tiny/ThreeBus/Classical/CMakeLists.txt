# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

gridkit_example_add_file(ThreeBusClassical.case.json)
gridkit_example_add_file(ThreeBusClassical.solver.json)
gridkit_example_add_file(ThreeBusClassical.ref.csv)

add_executable(ThreeBusClassical ThreeBusClassical.cpp)
target_link_libraries(
  ThreeBusClassical GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS ThreeBusClassical RUNTIME DESTINATION ${_install_path})

add_executable(ThreeBusClassicalJson ThreeBusClassicalJson.cpp)
target_link_libraries(
  ThreeBusClassicalJson GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS ThreeBusClassicalJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(ThreeBusClassical.json)

add_test(NAME ThreeBusClassical COMMAND ThreeBusClassical)
add_test(
  NAME ThreeBusClassical
  COMMAND DynamicSimulation ThreeBusClassical.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES ThreeBusClassical.case.json ThreeBusClassical.solver.json ThreeBusClassical.ref.csv
        DESTINATION ${_install_path})

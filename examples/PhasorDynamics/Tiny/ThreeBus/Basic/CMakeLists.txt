gridkit_example_current_install_path(_install_path)

add_executable(ThreeBusBasic ThreeBusBasic.cpp)
target_link_libraries(
  ThreeBusBasic GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS ThreeBusBasic RUNTIME DESTINATION ${_install_path})

add_executable(ThreeBusBasicJson ThreeBusBasicJson.cpp)
target_link_libraries(
  ThreeBusBasicJson
  GridKit::phasor_dynamics_systemmodel
  GridKit::solvers_dyn
  GridKit::utilities_cli_args
  GridKit::testing)
install(TARGETS ThreeBusBasicJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(ThreeBusBasic.case.json)
gridkit_example_add_file(ThreeBusBasic.ref.csv)

add_test(NAME ThreeBusBasic COMMAND ThreeBusBasic)
add_test(
  NAME ThreeBusBasicJson
  COMMAND ThreeBusBasicJson --case ThreeBusBasic.case.json --compare mon.csv ThreeBusBasic.ref.csv)

gridkit_example_add_file(ThreeBusBasic.solver.json)
add_test(
  NAME ThreeBusBasic_pdsim
  COMMAND DynamicSimulation ThreeBusBasic.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

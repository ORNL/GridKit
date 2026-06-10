gridkit_example_add_file(TwoBusTgov1.case.json)
gridkit_example_add_file(TwoBusTgov1.solver.json)
gridkit_example_add_file(TwoBusTgov1.ref.csv)

add_executable(TwoBusTgov1 TwoBusTgov1.cpp)
target_link_libraries(
  TwoBusTgov1 GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS TwoBusTgov1 RUNTIME DESTINATION ${_install_path})

add_executable(TwoBusTgov1Json TwoBusTgov1Json.cpp)
target_link_libraries(
  TwoBusTgov1Json GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS TwoBusTgov1Json RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(TwoBusTgov1.json)

add_test(NAME TwoBusTgov1 COMMAND TwoBusTgov1)
add_test(NAME TwoBusTgov1Json COMMAND TwoBusTgov1Json ${CMAKE_CURRENT_BINARY_DIR}/TwoBusTgov1.json)
add_test(
  NAME TwoBusTgov1Json
  COMMAND DynamicSimulation TwoBusTgov1.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES TwoBusTgov1.case.json TwoBusTgov1.solver.json TwoBusTgov1.ref.csv
        DESTINATION ${_install_path})

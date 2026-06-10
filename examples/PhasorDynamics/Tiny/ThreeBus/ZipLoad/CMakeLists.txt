gridkit_example_add_file(ThreeBusZipLoad.case.json)
gridkit_example_add_file(ThreeBusZipLoad.solver.json)
gridkit_example_add_file(ThreeBusZipLoad.ref.csv)

add_executable(ThreeBusZipLoadJson ThreeBusZipLoadJson.cpp)
target_link_libraries(
  ThreeBusZipLoadJson GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS ThreeBusZipLoadJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(ThreeBusZipLoad.json)


gridkit_example_current_install_path(_install_path)

add_executable(ThreeBusZipLoadJson ThreeBusZipLoadJson.cpp)
target_link_libraries(ThreeBusZipLoadJson
    GridKit::phasor_dynamics_components
    GridKit::solvers_dyn)
install(TARGETS ThreeBusZipLoadJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(ThreeBusZipLoad.json)

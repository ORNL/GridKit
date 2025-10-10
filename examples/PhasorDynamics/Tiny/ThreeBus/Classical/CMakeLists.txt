# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

gridkit_example_current_install_path(_install_path)

add_executable(ThreeBusClassical ThreeBusClassical.cpp)
target_link_libraries(ThreeBusClassical 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS ThreeBusClassical RUNTIME DESTINATION ${_install_path})

add_executable(ThreeBusClassicalJson ThreeBusClassicalJson.cpp)
target_link_libraries(ThreeBusClassicalJson
    GRIDKIT::phasor_dynamics_components
    GRIDKIT::solvers_dyn)
target_include_directories(ThreeBusClassicalJson
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/nlohmann-json/include)
target_include_directories(ThreeBusClassicalJson
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/magic-enum/include)
install(TARGETS ThreeBusClassicalJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(ThreeBusClassical.json)

add_test(NAME ThreeBusClassical COMMAND ThreeBusClassical)
add_test(NAME ThreeBusClassicalJson
    COMMAND
    ThreeBusClassicalJson ${CMAKE_CURRENT_BINARY_DIR}/ThreeBusClassical.json)
add_test(NAME ThreeBusClassicalJson_no_arg
    COMMAND ThreeBusClassicalJson
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

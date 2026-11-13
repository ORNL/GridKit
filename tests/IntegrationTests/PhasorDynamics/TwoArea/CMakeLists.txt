add_executable(test_pd_two_area_hierarchical runTwoAreaHierarchicalTests.cpp)
target_link_libraries(
  test_pd_two_area_hierarchical
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

configure_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/TwoArea/TwoArea.case.json TwoArea.case.json
               COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/TwoArea/TwoArea.state.json
               TwoArea.state.json COPYONLY)
configure_file(TwoAreaHierarchical.case.json TwoAreaHierarchical.case.json COPYONLY)

add_test(NAME PhasorDynamicsIntegrationTwoAreaHierarchical COMMAND test_pd_two_area_hierarchical)

install(TARGETS test_pd_two_area_hierarchical)
install(FILES ${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/TwoArea/TwoArea.case.json
              ${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/TwoArea/TwoArea.state.json
              TwoAreaHierarchical.case.json DESTINATION share/gridkit/tests/PhasorDynamics/TwoArea)

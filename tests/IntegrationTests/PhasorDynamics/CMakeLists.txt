add_executable(test_pd_component_integration runPDComponentIntegrationTests.cpp)
target_link_libraries(
  test_pd_component_integration
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::testing)

add_test(NAME PhasorDynamicsComponentIntegrationTest COMMAND test_pd_component_integration)

install(TARGETS test_pd_component_integration)

if(TARGET SUNDIALS::idas)
  add_executable(test_pd_integration runPDIntegrationTests.cpp)
  target_link_libraries(test_pd_integration PRIVATE GridKit::phasor_dynamics_systemmodel
                                                    GridKit::solvers_dyn GridKit::testing)

  add_test(
    NAME PhasorDynamicsIntegrationTest
    COMMAND test_pd_integration
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples/PhasorDynamics/Tiny)

  install(TARGETS test_pd_integration)
endif()

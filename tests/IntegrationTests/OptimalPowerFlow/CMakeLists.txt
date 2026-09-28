add_executable(test_opf_steady_start runSteadyStartTests.cpp)
target_link_libraries(
  test_opf_steady_start
  PRIVATE GridKit::phasor_dynamics_systemmodel
          GridKit::model_state_data
          GridKit::solvers_dyn
          GridKit::testing)
target_include_directories(
  test_opf_steady_start
  PRIVATE ${CMAKE_SOURCE_DIR}/application/PhasorDynamics)

install(TARGETS test_opf_steady_start)

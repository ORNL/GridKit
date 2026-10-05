add_executable(test_pd_integration runPDIntegrationTests.cpp)
target_link_libraries(
  test_pd_integration
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

add_test(
  NAME PhasorDynamicsIntegrationTest
  COMMAND test_pd_integration
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples/PhasorDynamics/Validation/Toy)

if(GRIDKIT_ENABLE_ENZYME AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_executable(test_pd_partitioned PartitionedIntegration.cpp)
  target_include_directories(test_pd_partitioned
                             PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
  target_link_libraries(test_pd_partitioned PRIVATE GridKit::phasor_dynamics_systemmodel
                                                    GridKit::solvers_dyn GridKit::testing)

  add_test(
    NAME PhasorDynamicsPartitionedIntegrationTest
    COMMAND
      test_pd_partitioned
      ${PROJECT_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusPartitionedA.case.json
      ${PROJECT_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusPartitionedB.case.json
      ${PROJECT_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusPartitioned.case.json)

  add_test(
    NAME PhasorDynamicsTenGenPartitionedIntegrationTest
    COMMAND
      test_pd_partitioned
      ${PROJECT_SOURCE_DIR}/cases/PhasorDynamics/Toy/TenGenPartitionedA.case.json
      ${PROJECT_SOURCE_DIR}/cases/PhasorDynamics/Toy/TenGenPartitionedB.case.json
      ${PROJECT_SOURCE_DIR}/cases/PhasorDynamics/Toy/TenGenPartitioned.case.json 10.0 0.0005)

  install(TARGETS test_pd_partitioned)
endif()

install(TARGETS test_pd_integration)

add_subdirectory(TenGen)

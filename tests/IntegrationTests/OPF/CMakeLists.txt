add_executable(test_opf_ipopt runIpoptSystemTests.cpp)
target_link_libraries(
  test_opf_ipopt
  PRIVATE GridKit::model_opf GridKit::solvers_opt_ipopt)
set_property(TARGET test_opf_ipopt PROPERTY LINKER_LANGUAGE Fortran)

add_test(
  NAME OPFIpoptExactTest
  COMMAND
    ${CMAKE_COMMAND} -DTEST_EXECUTABLE=$<TARGET_FILE:test_opf_ipopt> -P
    ${CMAKE_CURRENT_SOURCE_DIR}/CheckIpoptTest.cmake)

install(TARGETS test_opf_ipopt RUNTIME DESTINATION bin)

add_executable(test_opf_case_catalog runCaseCatalogTests.cpp)
target_link_libraries(
  test_opf_case_catalog
  PRIVATE GridKit::model_opf GridKit::solvers_opt_ipopt)
set_property(TARGET test_opf_case_catalog PROPERTY LINKER_LANGUAGE Fortran)

function(add_opf_catalog_case case_name state_subdirectory)
  add_test(
    NAME OPFCaseCatalog_${case_name}
    COMMAND
      test_opf_case_catalog ${PROJECT_SOURCE_DIR}/cases/OPF/${case_name}.opf.json
      ${PROJECT_SOURCE_DIR}/examples/PhasorDynamics/${state_subdirectory}/${case_name}.state.json
      solve)
endfunction()

add_opf_catalog_case(newengland Medium/NewEngland)
add_opf_catalog_case(hawaii Medium/Hawaii)
add_opf_catalog_case(illinois Large/Illinois)
add_opf_catalog_case(wecc Large/WECC)
add_opf_catalog_case(texas Large/Texas)

install(TARGETS test_opf_case_catalog RUNTIME DESTINATION bin)

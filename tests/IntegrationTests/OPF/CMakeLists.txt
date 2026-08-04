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

function(add_opf_catalog_case case_name)
  add_test(
    NAME OPFCaseCatalog_${case_name}
    COMMAND
      test_opf_case_catalog ${PROJECT_SOURCE_DIR}/cases/OPF/${case_name}.opf.json
      ${PROJECT_SOURCE_DIR}/cases/PhasorDynamics/${case_name}/${case_name}.state.json
      solve)
endfunction()

add_opf_catalog_case(IEEE39)
add_opf_catalog_case(Hawaii)
add_opf_catalog_case(ACTIVSg200)
add_opf_catalog_case(WECC240)
add_opf_catalog_case(ACTIVSg2000)

install(TARGETS test_opf_case_catalog RUNTIME DESTINATION bin)

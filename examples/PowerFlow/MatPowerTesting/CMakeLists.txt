
set(test_sources
  test_parse_gen_row
  test_parse_bus_row
  test_parse_branch_row
  test_parse_gencost_row
  test_parse_matpower
  )

foreach(src ${test_sources})
  add_executable(${src} ${src}.cpp)
  target_link_libraries(${src} PUBLIC GRIDKIT::Utilities)
  add_test(
    NAME ${src}
    COMMAND $<TARGET_FILE:${src}>
  )
endforeach()

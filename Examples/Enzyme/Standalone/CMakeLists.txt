enzyme_add_executable(
  NAME EnzymeStandaloneCheck
  SOURCES main.cpp
  )

add_test(NAME "EnzymeStandaloneCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneCheck)

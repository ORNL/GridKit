enzyme_add_executable(
  NAME EnzymeStandaloneScalarCheck
  SOURCES EnzymeScalar.cpp
  )

enzyme_add_executable(
  NAME EnzymeStandaloneVectorCheck
  SOURCES EnzymeVector.cpp
  )

add_test(NAME "EnzymeStandaloneScalarCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneScalarCheck)
add_test(NAME "EnzymeStandaloneVectorCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneVectorCheck)

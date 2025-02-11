enzyme_add_executable(
  NAME EnzymeLibScalarCheck 
  SOURCES EnzymeScalar.cpp ScalarModel.cpp
)

add_test(NAME "EnzymeLibScalarCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeLibScalarCheck)

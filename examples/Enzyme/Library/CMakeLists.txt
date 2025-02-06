enzyme_add_executable(NAME EnzymeLibCheck SOURCES main.cpp Library.cpp)

add_test(NAME "EnzymeLibCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeLibCheck)

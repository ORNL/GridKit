enzyme_add_executable(NAME EnzymeLibCheck SOURCES main.cpp library.cpp)
#install(TARGETS ${CMAKE_CURRENT_BINARY_DIR}/EnzymeLibCheck DESTINATION bin)

add_test(NAME "EnzymeLibCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeLibCheck)

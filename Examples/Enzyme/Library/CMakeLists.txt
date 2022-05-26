enzyme_add_executable(NAME libcheck SOURCES main.cpp library.cpp)
#install(TARGETS ${CMAKE_CURRENT_BINARY_DIR}/libcheck DESTINATION bin)

add_test(NAME "libcheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/libcheck)

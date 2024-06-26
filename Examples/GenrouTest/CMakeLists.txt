add_executable(genroutest GenrouTest.cpp)
target_link_libraries(genroutest GRIDKIT::bus GRIDKIT::genrou GRIDKIT::solvers_opt GRIDKIT::solvers_dyn)
install(TARGETS genroutest DESTINATION bin)

add_test(NAME GenrouTest COMMAND $<TARGET_FILE:genroutest>)
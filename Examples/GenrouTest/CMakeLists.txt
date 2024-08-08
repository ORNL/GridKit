add_executable(genroutest GenrouTest.cpp)
add_executable(genrougovtest GenrouGovTest.cpp)
target_link_libraries(genroutest GRIDKIT::bus GRIDKIT::genrou GRIDKIT::solvers_opt GRIDKIT::solvers_dyn)
target_link_libraries(genrougovtest GRIDKIT::bus GRIDKIT::genrou GRIDKIT::solvers_opt GRIDKIT::solvers_dyn GRIDKIT::governor GRIDKIT::load)
install(TARGETS genroutest DESTINATION bin)
install(TARGETS genrougovtest DESTINATION bin)

add_test(NAME GenrouTest COMMAND $<TARGET_FILE:genroutest>)
add_test(NAME GenrouGovTest COMMAND $<TARGET_FILE:genrougovtest>)
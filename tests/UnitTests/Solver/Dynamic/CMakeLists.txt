# 

add_executable(test_ida runIdaTests.cpp)
target_link_libraries(test_ida GRIDKIT::solvers_dyn)
                                

# add_test(NAME IDATest COMMAND $<TARGET_FILE:test_ida>)

install(TARGETS test_ida 
        RUNTIME DESTINATION bin)

# 

if(TARGET SUNDIALS::idas)
    add_executable(test_ida runIdaTests.cpp)
    target_link_libraries(test_ida GridKit::solvers_dyn GridKit::testing)
endif()

if (GRIDKIT_ENABLE_RESOLVE)
    add_executable(test_ros runRosenbrockTests.cpp)
    target_link_libraries(test_ros GridKit::rosenbrock GridKit::testing)
endif()
                                

# add_test(NAME IDATest COMMAND $<TARGET_FILE:test_ida>)

install(TARGETS test_ida 
        RUNTIME DESTINATION bin)

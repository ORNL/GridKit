# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]


add_executable(test_dependency_tracking runDependencyTrackingTests.cpp)
target_link_libraries(test_dependency_tracking GridKit::phasor_dynamics_bus)
                                

add_test(NAME DependencyTrackingTest COMMAND $<TARGET_FILE:test_dependency_tracking>)

install(TARGETS test_dependency_tracking 
        RUNTIME DESTINATION bin)

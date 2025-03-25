#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

# add_executable(grid3bus Grid3Bus.cpp)
# target_link_libraries(grid3bus GRIDKIT::minigrid GRIDKIT::solvers_steady)
# install(TARGETS grid3bus RUNTIME DESTINATION bin)

add_executable(grid3bus Grid3BusSys.cpp)
target_link_libraries(grid3bus GRIDKIT::minigrid 
                               GRIDKIT::bus
                               GRIDKIT::generator
                               GRIDKIT::branch
                               GRIDKIT::load
                               GRIDKIT::solvers_steady)
install(TARGETS grid3bus RUNTIME DESTINATION bin)

add_test(NAME Grid3Bus    COMMAND $<TARGET_FILE:grid3bus>)

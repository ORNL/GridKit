#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

# add_executable(grid3bus Grid3Bus.cpp)
# target_link_libraries(grid3bus GridKit::minigrid GridKit::solvers_steady)
# install(TARGETS grid3bus RUNTIME DESTINATION bin)

add_executable(grid3bus Grid3BusSys.cpp)
target_link_libraries(
  grid3bus
  GridKit::minigrid
  GridKit::bus
  GridKit::generator
  GridKit::branch
  GridKit::load
  GridKit::solvers_steady)
install(TARGETS grid3bus RUNTIME DESTINATION bin)

add_test(NAME Grid3Bus COMMAND $<TARGET_FILE:grid3bus>)

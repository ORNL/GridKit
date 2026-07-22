cmake_minimum_required(VERSION 3.20)
project(SolvePF)

find_package(GridKit REQUIRED)

add_executable(solve_pf solve_pf.cpp)
target_link_libraries(
  solve_pf
  bus
  generator
  branch
  load
  solvers_steady)
install(TARGETS solve_pf RUNTIME DESTINATION bin)

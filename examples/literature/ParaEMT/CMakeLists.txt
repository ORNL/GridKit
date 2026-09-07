cmake_minimum_required(VERSION 3.20)
project(ParaEMTGridKit LANGUAGES C CXX)
find_package(GridKit CONFIG REQUIRED)
get_filename_component(GRIDKIT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
add_executable(paraemt_9bus gridkit_9bus.cpp)
target_compile_features(paraemt_9bus PRIVATE cxx_std_20)
target_include_directories(
  paraemt_9bus
  PRIVATE "${GRIDKIT_SOURCE_DIR}" "${GRIDKIT_SOURCE_DIR}/third-party/nlohmann-json/include"
          "${GRIDKIT_SOURCE_DIR}/third-party/magic-enum/include")
target_link_libraries(
  paraemt_9bus
  PRIVATE emt_systemmodel solvers_dyn testing)

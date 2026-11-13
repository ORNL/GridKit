add_executable(TransientSimulation TransientSimulation.cpp)
target_link_libraries(
  TransientSimulation
  PUBLIC GridKit::emt_systemmodel GridKit::solvers_dyn GridKit::Utilities)
target_include_directories(
  TransientSimulation
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

install(
  TARGETS TransientSimulation
  EXPORT gridkit-targets
  RUNTIME)

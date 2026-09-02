# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#]]

add_executable(EMTDynamicSimulation DynamicSimulation.cpp)
target_link_libraries(
  EMTDynamicSimulation
  PUBLIC GridKit::emt_systemmodel
         GridKit::solvers_dyn
         GridKit::Utilities
         GridKit::testing)
target_include_directories(
  EMTDynamicSimulation
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
          ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

install(
  TARGETS EMTDynamicSimulation
  EXPORT gridkit-targets
  RUNTIME)

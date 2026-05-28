add_executable(DynamicSimulation DynamicSimulation.cpp)
target_link_libraries(DynamicSimulation
  PUBLIC
  GridKit::phasor_dynamics_components
  GridKit::solvers_dyn
  GridKit::Utilities
  GridKit::testing)
target_include_directories(DynamicSimulation PRIVATE
  ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
  ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

add_executable(ContingencyAnalysis ContingencyAnalysis.cpp)
target_link_libraries(ContingencyAnalysis
  PUBLIC
  GridKit::phasor_dynamics_components
  GridKit::solvers_dyn
  GridKit::Utilities
  GridKit::testing)
target_include_directories(ContingencyAnalysis PRIVATE
  ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
  ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

if(GRIDKIT_ENABLE_OPENMP)
  target_link_libraries(ContingencyAnalysis PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS DynamicSimulation EXPORT gridkit-targets RUNTIME)
install(TARGETS ContingencyAnalysis EXPORT gridkit-targets RUNTIME)

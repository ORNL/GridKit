add_executable(PDSim PDSim.cpp)
target_link_libraries(PDSim
  PUBLIC
  GridKit::phasor_dynamics_components
  GridKit::solvers_dyn
  GridKit::Utilities
  GridKit::testing)
target_include_directories(PDSim PRIVATE
  ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
  ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

add_executable(ContingencyAnalysis ContingencyAnalysis.cpp CTGStudy.cpp)
target_link_libraries(ContingencyAnalysis
  PUBLIC
  GridKit::phasor_dynamics_components
  GridKit::solvers_dyn
  GridKit::Utilities
  GridKit::testing)
target_include_directories(ContingencyAnalysis PRIVATE
  ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
  ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

option(GridKit_ENABLE_OPENMP "Use OpenMP in contingency analysis application" ON)
if(GridKit_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(ContingencyAnalysis PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS PDSim EXPORT gridkit-targets RUNTIME)
install(TARGETS ContingencyAnalysis EXPORT gridkit-targets RUNTIME)

add_executable(DynamicSimulation DynamicSimulation.cpp)
target_link_libraries(
  DynamicSimulation
  PUBLIC GridKit::phasor_dynamics_systemmodel
         GridKit::phasor_dynamics_partition
         GridKit::solvers_dyn
         GridKit::Utilities
         GridKit::testing)
target_include_directories(
  DynamicSimulation
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

add_executable(ContingencyAnalysis ContingencyAnalysis.cpp)
target_link_libraries(
  ContingencyAnalysis
  PUBLIC GridKit::phasor_dynamics_systemmodel
         GridKit::solvers_dyn
         GridKit::Utilities
         GridKit::testing)
target_include_directories(
  ContingencyAnalysis
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

if(GRIDKIT_ENABLE_OPENMP)
  target_link_libraries(ContingencyAnalysis PUBLIC OpenMP::OpenMP_CXX)
endif()

if(GRIDKIT_ENABLE_THREADS)
  target_link_libraries(ContingencyAnalysis PUBLIC Threads::Threads)
endif()

# Partitioned needs IDA and Ida needs IDAS, so they can't share an executable
if(GRIDKIT_ENABLE_SUNDIALS_PARTITIONED_DAE)
  add_executable(PartitionedDynamicSimulation PartitionedDynamicSimulation.cpp)
  target_link_libraries(
    PartitionedDynamicSimulation
    PUBLIC GridKit::phasor_dynamics_systemmodel
           GridKit::phasor_dynamics_partition
           GridKit::solvers_partitioned
           GridKit::Utilities
           GridKit::testing)
  # AnalysisUtilities.hpp includes Ida.hpp, which includes the KLU header
  # TODO: Move IdaConsistentICType to its own header and remove this link
  if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
    target_link_libraries(PartitionedDynamicSimulation PUBLIC SUNDIALS::sunlinsolklu)
  endif()
  target_include_directories(PartitionedDynamicSimulation
                             PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
  install(
    TARGETS PartitionedDynamicSimulation
    EXPORT gridkit-targets
    RUNTIME)
endif()

install(
  TARGETS DynamicSimulation
  EXPORT gridkit-targets
  RUNTIME)
install(
  TARGETS ContingencyAnalysis
  EXPORT gridkit-targets
  RUNTIME)

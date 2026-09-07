if(GRIDKIT_ENABLE_ENZYME)
  set(_source LineDistributedEnzyme.cpp)
else()
  set(_source LineDistributed.cpp)
endif()

gridkit_add_library(
  emt_line_distributed
  SOURCES ${_source}
  HEADERS LineDistributed.hpp LineDistributedData.hpp LineDistributedDataJSONParser.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_propagation)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_line_distributed)

gridkit_add_library(
  emt_line_distributed_dependency_tracking
  SOURCES LineDistributedDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_propagation_dependency_tracking)

target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_line_distributed_dependency_tracking)

if(GRIDKIT_ENABLE_ENZYME)
  target_link_libraries(emt_line_distributed PRIVATE ClangEnzymeFlags)
  target_compile_options(
    emt_line_distributed
    PRIVATE -mllvm
            -enzyme-auto-sparsity=1
            -fno-math-errno
            -fno-vectorize)
endif()

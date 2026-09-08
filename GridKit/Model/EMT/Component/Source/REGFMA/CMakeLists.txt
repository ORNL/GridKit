set(_install_headers Regfma.hpp RegfmaData.hpp RegfmaDataJSONParser.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    emt_regfma
    SOURCES RegfmaEnzyme.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::emt_core
      PUBLIC
      GridKit::emt_signal
      PRIVATE
      ClangEnzymeFlags
    COMPILE_OPTIONS
      PRIVATE
      -mllvm
      # Enzyme 16 stalls in automatic sparsity optimization for this residual.
      # SparseJacobian still lowers its explicit sparse load/store callbacks.
      -enzyme-auto-sparsity=0
      -fno-math-errno
      -fno-vectorize)
else()
  gridkit_add_library(
    emt_regfma
    SOURCES Regfma.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::emt_core
      PUBLIC
      GridKit::emt_signal)
endif()

gridkit_add_library(
  emt_regfma_dependency_tracking
  SOURCES RegfmaDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_regfma)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_regfma_dependency_tracking)

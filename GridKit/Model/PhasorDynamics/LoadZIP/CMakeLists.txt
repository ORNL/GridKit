
set(_install_headers
  LoadZIP.hpp
  LoadDataZIP.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_loadzip
    SOURCES
      LoadZIPEnzyme.cpp
    HEADERS
      ${_install_headers}
    INCLUDE_DIRECTORIES
      PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PRIVATE ClangEnzymeFlags
      PUBLIC GridKit::phasor_dynamics_core
    COMPILE_OPTIONS
      PRIVATE -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_loadzip
    SOURCES
      LoadZIP.cpp
    INCLUDE_DIRECTORIES
      PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC GridKit::phasor_dynamics_core
    HEADERS
      ${_install_headers})
endif()

gridkit_add_library(phasor_dynamics_loadzip_dependency_tracking
  SOURCES 
    LoadZIPDependencyTracking.cpp
  INCLUDE_DIRECTORIES 
    PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_loadzip)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_loadzip_dependency_tracking)

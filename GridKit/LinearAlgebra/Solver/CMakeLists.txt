if(GRIDKIT_ENABLE_RESOLVE)
  gridkit_add_library(
    resolve_system_solver
    SOURCES ResolveSystemSolver.cpp
    HEADERS LinearSolver.hpp ResolveSystemSolver.hpp
    LINK_LIBRARIES PRIVATE ReSolve::ReSolve)
endif()

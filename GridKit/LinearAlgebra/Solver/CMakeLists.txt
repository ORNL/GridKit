install(FILES LinearSolver.hpp DESTINATION include/GridKit/LinearAlgebra/Solver)

if(GRIDKIT_ENABLE_RESOLVE)
  gridkit_add_library(
    resolve_system_solver
    SOURCES ResolveSystemSolver.cpp
    HEADERS ResolveSystemSolver.hpp
    LINK_LIBRARIES PRIVATE ReSolve::ReSolve)
endif()

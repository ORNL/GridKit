gridkit_add_library(
  solvers_opt_ipopt
  SOURCES IpoptAdapter.cpp
  HEADERS IpoptAdapter.hpp
  LINK_LIBRARIES
    PUBLIC
    IPOPT
    PUBLIC
    GridKit::optimization_core)

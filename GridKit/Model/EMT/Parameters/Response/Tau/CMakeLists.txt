gridkit_add_library(
  emt_tau
  SOURCES Tau.cpp
  HEADERS Tau.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_gamma
    GridKit::emt_path)

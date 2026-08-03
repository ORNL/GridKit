gridkit_add_library(
  emt_h
  SOURCES H.cpp
  HEADERS H.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_gamma
    GridKit::emt_path)

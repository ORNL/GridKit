gridkit_add_library(
  emt_carson
  SOURCES Carson.cpp
  HEADERS Carson.hpp CarsonData.hpp
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_tower)

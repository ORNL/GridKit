gridkit_add_library(
  emt_skin_effect
  SOURCES SkinEffect.cpp
  HEADERS SkinEffect.hpp
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_conductor)

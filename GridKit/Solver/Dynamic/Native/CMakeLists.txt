if(GRIDKIT_ENABLE_RESOLVE)
  gridkit_add_library(
    inf_norm
    SOURCES InfNorm.cpp
    HEADERS InfNorm.hpp
    LINK_LIBRARIES
      PRIVATE
      ReSolve::ReSolve
      PUBLIC
      GridKit::dense_vector)
endif()

gridkit_add_library(
  adaptive_step
  SOURCES AdaptiveStep.cpp
  HEADERS AdaptiveStep.hpp)

gridkit_add_library(
  fixed_step
  SOURCES FixedStep.cpp
  HEADERS FixedStep.hpp)

gridkit_add_library(
  utilities_logger
  SOURCES Logger.cpp
  HEADERS Logger.hpp
  LINK_LIBRARIES PRIVATE GridKit::definitions)

target_link_libraries(
  Utilities
  INTERFACE GridKit::utilities_logger)

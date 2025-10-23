gridkit_add_library(utilities_logger
  SOURCES
    Logger.cpp
  OUTPUT_NAME
    gridkit_utilities_logger)

target_link_libraries(Utilities INTERFACE GRIDKIT::utilities_logger)

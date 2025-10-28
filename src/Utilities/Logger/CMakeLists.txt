gridkit_add_library(utilities_logger
  SOURCES
    Logger.cpp)

target_link_libraries(Utilities INTERFACE GRIDKIT::utilities_logger)

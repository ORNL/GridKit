gridkit_add_library(utilities_cli_options
  SOURCES
    CliOptions.cpp
  HEADERS
    CliOptions.hpp
  LINK_LIBRARIES
    PUBLIC GridKit::utilities_logger
  )

target_link_libraries(Utilities INTERFACE GridKit::utilities_cli_options)

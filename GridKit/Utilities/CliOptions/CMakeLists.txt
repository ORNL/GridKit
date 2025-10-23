gridkit_add_library(utilities_cli_options
  SOURCES
    CliOptions.cpp
  HEADERS
    CliOptions.hpp)

target_link_libraries(Utilities INTERFACE GridKit::utilities_cli_options)

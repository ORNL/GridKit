gridkit_add_library(utilities_cli_options
  SOURCES
    CliOptions.cpp
  OUTPUT_NAME
    gridkit_utilities_cli_options)

target_link_libraries(Utilities INTERFACE GRIDKIT::utilities_cli_options)

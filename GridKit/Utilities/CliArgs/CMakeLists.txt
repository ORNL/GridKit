gridkit_add_library(utilities_cli_args
  SOURCES
    CliArgs.cpp
  HEADERS
    CliArgs.hpp
    ArgValue.hpp
    ArgVector.hpp
    Option.hpp
  LINK_LIBRARIES
    PUBLIC GridKit::utilities_logger
  )

target_link_libraries(Utilities INTERFACE GridKit::utilities_cli_args)

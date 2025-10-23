add_library(Utilities INTERFACE)
target_include_directories(Utilities INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)

add_library(GRIDKIT::Utilities ALIAS Utilities)

add_subdirectory(Logger)
add_subdirectory(CliOptions)

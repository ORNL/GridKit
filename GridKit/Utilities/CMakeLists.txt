add_library(Utilities INTERFACE)
target_include_directories(Utilities INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)

add_library(GridKit::Utilities ALIAS Utilities)

add_subdirectory(Logger)
add_subdirectory(CliOptions)

install(
  FILES
    Colors.hpp
    FileIO.hpp
    MapFromCOO.hpp
  DESTINATION
    include/GridKit/Utilities)

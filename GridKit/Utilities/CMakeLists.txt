add_library(Utilities INTERFACE)
target_include_directories(
  Utilities
  INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_compile_features(Utilities INTERFACE cxx_std_20)
target_link_libraries(
  Utilities
  INTERFACE GridKit::magic_enum)

add_library(GridKit::Utilities ALIAS Utilities)

add_subdirectory(Logger)
add_subdirectory(CliArgs)

install(TARGETS Utilities EXPORT gridkit-targets)
install(
  FILES Colors.hpp
        Enum.hpp
        Errors.hpp
        FileIO.hpp
        MapFromCsr.hpp
        Stream.hpp
        String.hpp
  DESTINATION include/GridKit/Utilities)

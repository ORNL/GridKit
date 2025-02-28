add_library(Utilities INTERFACE)
target_include_directories(Utilities INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(GRIDKIT::Utilities ALIAS Utilities)


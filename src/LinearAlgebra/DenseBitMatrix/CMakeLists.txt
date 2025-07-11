add_library(DenseBitMatrix INTERFACE)
include_directories(DenseBitMatrix INTERFACE ${CMAKE_CURRENT_LIST_DIR})

add_library(GRIDKIT::DenseBitMatrix ALIAS DenseBitMatrix)

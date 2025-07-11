add_library(CsrMatrix INTERFACE)
include_directories(CsrMatrix INTERFACE ${CMAKE_CURRENT_LIST_DIR})

add_library(GRIDKIT::CsrMatrix ALIAS CsrMatrix)

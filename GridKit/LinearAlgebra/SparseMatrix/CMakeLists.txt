gridkit_add_library(
  sparse_matrix
  SOURCES CooMatrix.cpp CsrMatrix.cpp
  HEADERS CooMatrix.hpp CsrMatrix.hpp
  LINK_LIBRARIES GridKit::utilities_logger GridKit::cpu_backend)

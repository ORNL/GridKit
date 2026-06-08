gridkit_add_library(
  dense_vector
  SOURCES Vector.cpp
  HEADERS Vector.hpp VectorView.hpp
  LINK_LIBRARIES GridKit::utilities_logger)

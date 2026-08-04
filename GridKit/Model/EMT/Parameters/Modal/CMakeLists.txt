# Identity-preserving modal observation of the propagation matrix.

gridkit_add_library(
  emt_modal
  SOURCES ModalDecomposition.cpp
  HEADERS ModalDecomposition.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PRIVATE
    Eigen3::Eigen)

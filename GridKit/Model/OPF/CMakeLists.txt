gridkit_add_library(
  model_opf
  SOURCES System.cpp
  HEADERS Branch.hpp
          Bus.hpp
          Component.hpp
          ComponentConstraints.hpp
          ComponentVariables.hpp
          Generator.hpp
          Load.hpp
          Shunt.hpp
          System.hpp
          SystemData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::model_state_data
    PUBLIC
    GridKit::dense_vector
    PUBLIC
    GridKit::sparse_matrix)

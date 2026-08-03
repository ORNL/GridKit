# [[
#  EMT subsystem foundation: shared wiring (Signal), model contract (Element),
#  and assembly utilities. Requires Eigen.
# ]]

gridkit_add_library(
  emt_core
  INTERFACE_TARGET
  HEADERS Constants.hpp
          Signal.hpp
          Element.hpp
          SystemAssembly.hpp
  LINK_LIBRARIES
    INTERFACE
    GridKit::definitions
    INTERFACE
    GridKit::utilities_logger
    INTERFACE
    GridKit::sparse_matrix)

add_subdirectory(Parameters)

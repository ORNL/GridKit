# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  minigrid
  SOURCES MiniGrid.cpp
  HEADERS MiniGrid.hpp
  LINK_LIBRARIES GridKit::dense_vector)

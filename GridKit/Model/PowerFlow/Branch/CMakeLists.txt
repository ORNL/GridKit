# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  branch
  SOURCES Branch.cpp
  HEADERS Branch.hpp
  LINK_LIBRARIES GridKit::dense_vector)

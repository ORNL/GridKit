# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  load
  SOURCES Load.cpp
  HEADERS Load.hpp
  LINK_LIBRARIES GridKit::dense_vector)

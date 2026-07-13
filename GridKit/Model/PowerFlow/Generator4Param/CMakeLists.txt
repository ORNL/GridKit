#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  generator4param
  SOURCES Generator4Param.cpp
  HEADERS Generator4Param.hpp
  LINK_LIBRARIES GridKit::dense_vector)

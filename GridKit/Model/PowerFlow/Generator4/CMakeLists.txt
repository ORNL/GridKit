#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  generator4
  SOURCES Generator4.cpp
  HEADERS Generator4.hpp
  LINK_LIBRARIES GridKit::dense_vector)

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(generator
  SOURCES
    GeneratorSlack.cpp
    GeneratorPV.cpp
    GeneratorPQ.cpp
  OUTPUT_NAME
    gridkit_generator)

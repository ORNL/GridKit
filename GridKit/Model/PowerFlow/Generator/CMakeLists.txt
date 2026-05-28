# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  generator
  SOURCES GeneratorSlack.cpp GeneratorPV.cpp GeneratorPQ.cpp
  HEADERS GeneratorBase.hpp
          GeneratorFactory.hpp
          GeneratorPQ.hpp
          GeneratorPV.hpp
          GeneratorSlack.hpp)

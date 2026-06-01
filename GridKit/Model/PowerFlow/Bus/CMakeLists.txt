#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  bus
  SOURCES BusPQ.cpp BusPV.cpp BusSlack.cpp
  HEADERS BaseBus.hpp
          BusFactory.hpp
          BusPQ.hpp
          BusPV.hpp
          BusSlack.hpp)

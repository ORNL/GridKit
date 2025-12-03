#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_subdirectory(PhasorDynamics)
add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)

install(FILES
  Evaluator.hpp
  VariableMonitorBase.hpp
  VariableMonitor.hpp
  DESTINATION include/GridKit/Model)

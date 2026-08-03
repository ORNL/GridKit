#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_subdirectory(PhasorDynamics)
add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)

if(GRIDKIT_ENABLE_EIGEN)
  add_subdirectory(EMT)
endif()

install(
  FILES Evaluator.hpp
        LogEvaluator.hpp
        VariableMonitor.hpp
        VariableMonitorController.hpp
  DESTINATION include/GridKit/Model)

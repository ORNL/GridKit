#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_subdirectory(PhasorDynamics)
add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(OptimalPowerFlow)
endif()

install(
  FILES ConfigurationChecks.hpp
        Evaluator.hpp
        OptimizationEvaluator.hpp
        ParameterReader.hpp
        VariableMonitor.hpp
        VariableMonitorController.hpp
  DESTINATION include/GridKit/Model)

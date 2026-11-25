#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_subdirectory(EMT)
add_subdirectory(PhasorDynamics)
add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)
add_subdirectory(OptimalPowerFlow)

install(
  FILES ConfigurationChecks.hpp
        DataReader.hpp
        Evaluator.hpp
        OptimizationEvaluator.hpp
        VariableMonitor.hpp
        VariableMonitorController.hpp
  DESTINATION include/GridKit/Model)

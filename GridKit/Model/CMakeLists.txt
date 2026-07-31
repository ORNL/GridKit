#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  model_state_data
  SOURCES StateData.cpp
  HEADERS StateData.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

add_subdirectory(PhasorDynamics)
add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)
add_subdirectory(OPF)

install(FILES Evaluator.hpp VariableMonitor.hpp VariableMonitorController.hpp
        DESTINATION include/GridKit/Model)

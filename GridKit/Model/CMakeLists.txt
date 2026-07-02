#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  model_variable_monitor
  SOURCES VariableMonitorArrowOutput.cpp
  HEADERS VariableMonitorArrowOutput.hpp
  LINK_LIBRARIES PRIVATE GridKit::definitions)
if(GRIDKIT_ENABLE_ARROW)
  target_link_libraries(model_variable_monitor PRIVATE ${GRIDKIT_ARROW_TARGET})
endif()

add_subdirectory(PhasorDynamics)
add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)

install(FILES Evaluator.hpp VariableMonitor.hpp VariableMonitorController.hpp
        DESTINATION include/GridKit/Model)

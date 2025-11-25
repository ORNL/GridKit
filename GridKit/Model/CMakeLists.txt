#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

add_subdirectory(PhasorDynamics)
add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)

install(FILES Evaluator.hpp DESTINATION include/GridKit/Model)
install(FILES EvaluatorMixins.hpp DESTINATION include/GridKit/Model)

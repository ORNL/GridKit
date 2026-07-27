# CMake generated Testfile for 
# Source directory: /home/isatkaus/gridkit/examples/PhasorDynamics/Medium/Hawaii
# Build directory: /home/isatkaus/gridkit/build/examples/PhasorDynamics/Medium/Hawaii
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(hawaii_pdsim "/home/isatkaus/gridkit/build/application/PhasorDynamics/DynamicSimulation" "hawaii.solver.json")
set_tests_properties(hawaii_pdsim PROPERTIES  WORKING_DIRECTORY "/home/isatkaus/gridkit/build/examples/PhasorDynamics/Medium/Hawaii" _BACKTRACE_TRIPLES "/home/isatkaus/gridkit/examples/PhasorDynamics/Medium/Hawaii/CMakeLists.txt;4;add_test;/home/isatkaus/gridkit/examples/PhasorDynamics/Medium/Hawaii/CMakeLists.txt;0;")

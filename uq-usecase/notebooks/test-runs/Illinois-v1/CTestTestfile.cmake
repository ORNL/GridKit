# CMake generated Testfile for 
# Source directory: /home/isatkaus/gridkit/examples/PhasorDynamics/Large/Illinois
# Build directory: /home/isatkaus/gridkit/build/examples/PhasorDynamics/Large/Illinois
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(illinois_pdsim "/home/isatkaus/gridkit/build/application/PhasorDynamics/DynamicSimulation" "illinois.solver.json")
set_tests_properties(illinois_pdsim PROPERTIES  WORKING_DIRECTORY "/home/isatkaus/gridkit/build/examples/PhasorDynamics/Large/Illinois" _BACKTRACE_TRIPLES "/home/isatkaus/gridkit/examples/PhasorDynamics/Large/Illinois/CMakeLists.txt;4;add_test;/home/isatkaus/gridkit/examples/PhasorDynamics/Large/Illinois/CMakeLists.txt;0;")

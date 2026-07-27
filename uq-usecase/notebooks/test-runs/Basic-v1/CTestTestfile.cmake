# CMake generated Testfile for 
# Source directory: /home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic
# Build directory: /home/isatkaus/gridkit/build/examples/PhasorDynamics/Tiny/ThreeBus/Basic
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ThreeBusBasic "/home/isatkaus/gridkit/build/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic")
set_tests_properties(ThreeBusBasic PROPERTIES  _BACKTRACE_TRIPLES "/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/CMakeLists.txt;19;add_test;/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/CMakeLists.txt;0;")
add_test(ThreeBusBasicJson "/home/isatkaus/gridkit/build/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson" "--case" "ThreeBusBasic.case.json" "--compare" "mon.csv" "ThreeBusBasic.ref.csv")
set_tests_properties(ThreeBusBasicJson PROPERTIES  _BACKTRACE_TRIPLES "/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/CMakeLists.txt;20;add_test;/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/CMakeLists.txt;0;")
add_test(ThreeBusBasic_pdsim "/home/isatkaus/gridkit/build/application/PhasorDynamics/DynamicSimulation" "ThreeBusBasic.solver.json")
set_tests_properties(ThreeBusBasic_pdsim PROPERTIES  WORKING_DIRECTORY "/home/isatkaus/gridkit/build/examples/PhasorDynamics/Tiny/ThreeBus/Basic" _BACKTRACE_TRIPLES "/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/CMakeLists.txt;25;add_test;/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/CMakeLists.txt;0;")

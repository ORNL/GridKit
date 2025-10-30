# Path where the consumer test code will be installed
set(_rel_consumer_path ${GRIDKIT_EXAMPLES_INSTALL_ROOT}/GridKitConsumer)
set(_abs_consumer_path ${CMAKE_INSTALL_PREFIX}/${_rel_consumer_path})

# Install script as example of how to build and test consumer project
install(PROGRAMS test.sh DESTINATION ${_rel_consumer_path})

# Install consumer code and set link libraries
set(_examples_dir ${CMAKE_SOURCE_DIR}/examples)
if(GRIDKIT_ENABLE_SUNDIALS)
  install(
    FILES ${_examples_dir}/PhasorDynamics/Tiny/TwoBus/Basic/TwoBusBasicJson.cpp
    DESTINATION ${_rel_consumer_path}
    RENAME consumer.cpp)
  install(
    FILES
      ${_examples_dir}/PhasorDynamics/Tiny/TwoBus/Basic/TwoBusBasic.hpp
      ${_examples_dir}/PhasorDynamics/Tiny/TwoBus/Basic/TwoBusBasic.json
    DESTINATION ${_rel_consumer_path})
  list(APPEND GRIDKIT_CONSUMER_TEST_LINK_LIBRARIES
    GRIDKIT::phasor_dynamics_components
    GRIDKIT::solvers_dyn)
else()
  install(
    FILES ${_examples_dir}/PowerElectronics/DistributedGeneratorTest/DGTest.cpp
    DESTINATION ${_rel_consumer_path}
    RENAME consumer.cpp)
  list(APPEND GRIDKIT_CONSUMER_TEST_LINK_LIBRARIES
    GRIDKIT::power_elec_disgen
    GRIDKIT::power_elec_microline
    GRIDKIT::power_elec_microload)
endif()

# Install CMakeLists.txt for example project
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GridKitConsumer/CMakeLists.txt.in
    ${CMAKE_CURRENT_BINARY_DIR}/GridKitConsumer/CMakeLists.txt
    @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GridKitConsumer/CMakeLists.txt
  DESTINATION ${_rel_consumer_path})

# Add target to build and test the example project after install
add_custom_target(test_install COMMAND
    ${_abs_consumer_path}/test.sh ${CMAKE_INSTALL_PREFIX} ${CMAKE_CXX_COMPILER})

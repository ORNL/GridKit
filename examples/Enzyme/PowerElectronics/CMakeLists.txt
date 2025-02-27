enzyme_add_executable(
  NAME EnzymePowerElectronicsCheck 
  SOURCES main.cpp
  LINK_LIBRARIES GRIDKIT::DenseMatrix GRIDKIT::power_elec_disgen
)

add_test(NAME "EnzymePowerElectronicsCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymePowerElectronicsCheck)

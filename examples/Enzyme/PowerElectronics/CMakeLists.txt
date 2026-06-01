add_executable(EnzymePowerElectronicsCheck main.cpp)
target_compile_options(EnzymePowerElectronicsCheck PUBLIC -fno-vectorize -ffast-math
                                                          -fno-unroll-loops)
target_link_libraries(
  EnzymePowerElectronicsCheck
  ClangEnzymeFlags
  GridKit::DenseMatrix
  GridKit::power_elec_disgen
  GridKit::Utilities)

add_test(NAME "EnzymePowerElectronicsCheck"
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymePowerElectronicsCheck)

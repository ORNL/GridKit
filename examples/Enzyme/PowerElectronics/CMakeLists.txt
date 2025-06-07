add_executable(EnzymePowerElectronicsCheck main.cpp)
target_compile_options(EnzymePowerElectronicsCheck PUBLIC -O2 -fno-vectorize -ffast-math -fno-unroll-loops)
target_link_libraries(EnzymePowerElectronicsCheck ClangEnzymeFlags GRIDKIT::DenseMatrix GRIDKIT::power_elec_disgen)

add_test(NAME "EnzymePowerElectronicsCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymePowerElectronicsCheck)

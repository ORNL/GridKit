if(NOT GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  return()
endif()

# This standalone EMT example reuses the existing signal-delay assembly.
add_executable(line_distributed EXCLUDE_FROM_ALL main.cpp LineCircuit.cpp)
add_executable(line_distributed_tests EXCLUDE_FROM_ALL tests.cpp LineCircuit.cpp)
foreach(target line_distributed line_distributed_tests)
  target_link_libraries(${target} PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
  target_include_directories(${target} PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
endforeach()

foreach(
  file
  demo.json
  physical.json
  bergeron.fit.json
  frequency-dependent.fit.json
  line_constants.py
  fit_line.py
  validate_physical.py
  plot_physical.py
  validate.py
  plot.py)
  gridkit_example_add_file(${file})
endforeach()

add_subdirectory(data/overhead)

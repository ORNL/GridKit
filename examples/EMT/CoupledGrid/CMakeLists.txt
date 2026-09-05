# Keep the parameter data and reproduction scripts together with the case.
foreach(
  _file
  CoupledGrid.case.json
  CoupledGrid.state.json
  low.solver.json
  middle.solver.json
  high.solver.json
  README.md
  build_case.py
  run.py
  analyze.py
  report.py
  pwm_analysis.py
  validate_case.py
  draw_oneline.py
  plot_results.py
  lines/README.md
  lines/model.json
  lines/fit_validation.json
  lines/openline.csv
  lines/response.csv
  lines/openline_sweep.rs
  lines/Cargo.lock
  lines/generate.py)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/CoupledGrid/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/CoupledGrid/${_file} COPYONLY)
  get_filename_component(_directory ${_file} DIRECTORY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/CoupledGrid/${_file}
          DESTINATION share/gridkit/cases/EMT/CoupledGrid/${_directory})
endforeach()

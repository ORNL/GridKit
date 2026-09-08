foreach(
  _file
  GFL.solver.json
  GFM.solver.json
  README.md
  run.py
  plot.py)
  configure_file(${_file} ${_file} COPYONLY)
  install(FILES ${_file} DESTINATION share/gridkit/examples/EMT/CurrentControl)
endforeach()

foreach(
  _file
  GFL.case.json
  GFM.case.json
  GFL.state.json
  GFM.state.json
  README.md)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/CurrentControl/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/CurrentControl/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/CurrentControl/${_file}
          DESTINATION share/gridkit/cases/EMT/CurrentControl)
endforeach()

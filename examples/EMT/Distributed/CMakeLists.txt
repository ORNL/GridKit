foreach(
  _file
  README.md
  requirements.txt
  fit.py
  run.py
  plot.py)
  gridkit_example_add_file(${_file})
endforeach()

foreach(_file README.md overhead.line.json north-american.catalog.json)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/Distributed/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/Distributed/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/Distributed/${_file}
          DESTINATION share/gridkit/cases/EMT/Distributed)
endforeach()

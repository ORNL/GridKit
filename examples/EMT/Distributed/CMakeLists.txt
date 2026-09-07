foreach(
  _file
  README.md
  requirements.txt
  fit.py
  run.py
  pulse.py
  pulse_plot.py
  pulse_geometry.tex
  pulse_circuits.tex
  plot.py)
  gridkit_example_add_file(${_file})
endforeach()

gridkit_example_add_file(${CMAKE_SOURCE_DIR}/docs/Figures/EMT/diagram-style.tex)

foreach(_file README.md overhead.line.json north-american.catalog.json)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/Distributed/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/Distributed/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/Distributed/${_file}
          DESTINATION share/gridkit/cases/EMT/Distributed)
endforeach()

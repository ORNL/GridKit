foreach(
  _file
  FaultClearing.solver.json
  run.py
  plot.py
  README.md)
  gridkit_example_add_file(${_file})
endforeach()

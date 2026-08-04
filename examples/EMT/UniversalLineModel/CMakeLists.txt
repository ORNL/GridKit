# Universal line model example: fit the committed line description and
# regenerate the model files in the build tree. Every mode is unwound
# by its own delay, so no residual winding constrains the band; the
# reduced band and loose target just keep the example quick.

gridkit_example_current_install_path(_install_path)

gridkit_example_add_file(unsymmetric-flat.line.json)

add_test(
  NAME UniversalLineModelExample
  COMMAND UniversalLineModel unsymmetric-flat.line.json -o . --h-target 0.2 --fmax 1e6)

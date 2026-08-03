# Universal line model example: fit the committed line description and
# regenerate the model files in the build tree. The propagation factor
# target reflects what factor-based fitting supports for a line with
# nearly degenerate aerial modes; see the application header note.

gridkit_example_current_install_path(_install_path)

gridkit_example_add_file(unsymmetric-flat.line.json)

add_test(
  NAME UniversalLineModelExample
  COMMAND UniversalLineModel unsymmetric-flat.line.json -o . --h-target 0.05)

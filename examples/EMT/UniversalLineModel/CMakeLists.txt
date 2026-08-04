# Universal line model example: fit the committed line description and
# regenerate the model files in the build tree. The propagation factor
# target reflects what factor-based fitting supports for a line with
# nearly degenerate aerial modes now that the sweep is error-controlled;
# see the application header note.

gridkit_example_current_install_path(_install_path)

gridkit_example_add_file(unsymmetric-flat.line.json)

add_test(
  NAME UniversalLineModelExample
  COMMAND UniversalLineModel unsymmetric-flat.line.json -o . --h-domain modal --h-target 0.2)

# The phase treatment carries the inter-mode delay spread as residual
# winding that grows with frequency, so its example runs on the band
# where one shared delay is apt; see the application README.
add_test(
  NAME UniversalLineModelPhaseExample
  COMMAND
    UniversalLineModel unsymmetric-flat.line.json -o phase --h-domain phase --h-target 0.2 --fmax
    1e6)

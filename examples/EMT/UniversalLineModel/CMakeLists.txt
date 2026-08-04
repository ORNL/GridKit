# Universal line model example: fit the committed line description and
# regenerate the model files in the build tree. One shared delay leaves
# the inter-mode spread as residual winding that grows with frequency,
# so the example runs on the band where a single delay is apt; see the
# application README.

gridkit_example_current_install_path(_install_path)

gridkit_example_add_file(unsymmetric-flat.line.json)

add_test(
  NAME UniversalLineModelExample
  COMMAND UniversalLineModel unsymmetric-flat.line.json -o . --h-target 0.2 --fmax 1e6)

# Preserve the legacy 39-fault study; the reusable IEEE39 case has one fault.
configure_file(newengland.case.json newengland.case.json COPYONLY)
configure_file(newengland.solver.json newengland.solver.json COPYONLY)

add_test(
  NAME newengland_ca
  COMMAND ContingencyAnalysis newengland.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

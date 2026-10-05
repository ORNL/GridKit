gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/ACTIVSg10k/ACTIVSg10k.case.json)
foreach(regions 4 8 16)
  gridkit_example_add_file(
    ${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/ACTIVSg10k/ACTIVSg10k.${regions}.partition.json)
endforeach()
gridkit_example_add_file(ACTIVSg10kPartitioned.solver.json)
gridkit_example_add_file(README.md)

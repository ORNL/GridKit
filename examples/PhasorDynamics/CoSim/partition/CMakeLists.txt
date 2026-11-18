gridkit_example_add_file(IEEE39.A.case.json)
gridkit_example_add_file(IEEE39.B.case.json)

foreach(region IN LISTS wecc_regions)
  gridkit_example_add_file(WECC240.${region}.case.json)
endforeach()

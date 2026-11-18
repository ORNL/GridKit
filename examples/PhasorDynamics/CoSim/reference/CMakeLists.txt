gridkit_example_add_file(IEEE39.A.ref.csv)
gridkit_example_add_file(IEEE39.B.ref.csv)

foreach(region IN LISTS wecc_regions)
  gridkit_example_add_file(WECC240.${region}.ref.csv)
endforeach()

add_executable(scalemicrogrid ScaleMicrogrid.cpp)
add_executable(scalemicrogridarbitrary ScaleMicrogridArbitrary.cpp)

target_link_libraries(
  scalemicrogrid GridKit::solvers_dyn power_elec_microgrid_network)

target_link_libraries(
  scalemicrogridarbitrary GridKit::solvers_dyn power_elec_microgrid_network)

add_test(NAME ScaleMicrogrid COMMAND $<TARGET_FILE:scalemicrogrid>)
add_test(NAME scalemicrogridarbitrary COMMAND $<TARGET_FILE:scalemicrogridarbitrary>)
install(TARGETS scalemicrogrid RUNTIME DESTINATION bin)
install(TARGETS scalemicrogridarbitrary RUNTIME DESTINATION bin)

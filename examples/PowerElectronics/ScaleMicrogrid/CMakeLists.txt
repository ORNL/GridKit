add_executable(scalemicrogrid ScaleMicrogrid.cpp)
add_executable(scalemicrogridarbitrary ScaleMicrogridArbitrary.cpp)

target_include_directories(
  scalemicrogrid
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_include_directories(
  scalemicrogridarbitrary
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(
  scalemicrogrid
  GridKit::power_elec_disgen
  GridKit::power_elec_microline
  GridKit::power_elec_microload
  GridKit::solvers_dyn
  GridKit::power_elec_microbusdq)

target_link_libraries(
  scalemicrogridarbitrary
  GridKit::power_elec_disgen
  GridKit::power_elec_microline
  GridKit::power_elec_microload
  GridKit::solvers_dyn
  GridKit::power_elec_microbusdq)

add_test(NAME ScaleMicrogrid COMMAND $<TARGET_FILE:scalemicrogrid>)
install(TARGETS scalemicrogrid RUNTIME DESTINATION bin)

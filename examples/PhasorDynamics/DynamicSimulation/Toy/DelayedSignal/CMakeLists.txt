gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/DelayedSignal.case.json)
gridkit_example_add_file(DelayedSignal.solver.json)

# Build only when reproducing the diagnostic gallery.
add_executable(signal_delay_gallery EXCLUDE_FROM_ALL gallery/gallery.cpp)
target_link_libraries(
  signal_delay_gallery
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
target_include_directories(
  signal_delay_gallery
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

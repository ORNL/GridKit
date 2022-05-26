enzyme_add_executable(
  NAME StandaloneEnzymeCheck
  SOURCES main.cpp
  )

  add_test(NAME "standalonecheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/StandaloneEnzymeCheck)

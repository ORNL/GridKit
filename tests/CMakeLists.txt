if(NOT BUILD_TESTING)
  return()
endif()

add_subdirectory(UnitTests)
add_subdirectory(IntegrationTests)

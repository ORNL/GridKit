if(NOT BUILD_TESTING)
    return()
endif()

add_subdirectory(UnitTests)

if(GRIDKIT_ENABLE_BENCHMARKS)
add_subdirectory(Benchmarks)
endif()
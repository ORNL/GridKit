#[[

@brief Build GridKit memory management utilities

@author Slaven Peles <peless@ornl.gov>

]]

add_subdirectory(cpu)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/MemoryUtils.hpp ${CMAKE_CURRENT_SOURCE_DIR}/MemoryUtils.tpp
        DESTINATION include/GridKit/MemoryUtilities)

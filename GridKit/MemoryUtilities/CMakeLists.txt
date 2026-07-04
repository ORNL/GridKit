#[[

@brief Build GridKit memory management utilities

@author Slaven Peles <peless@ornl.gov>

]]

add_subdirectory(cpu)
# add_subdirectory(cuda)
# add_subdirectory(hip)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/MemoryUtils.hpp ${CMAKE_CURRENT_SOURCE_DIR}/MemoryUtils.tpp
        DESTINATION include/GridKit/MemoryUtilities)

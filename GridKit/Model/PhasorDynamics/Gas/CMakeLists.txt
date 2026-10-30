# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

install(FILES GasBusBase.hpp DESTINATION include/GridKit/Model/PhasorDynamics/Gas)

add_subdirectory(GasBranch)
add_subdirectory(GasBus)
add_subdirectory(GasBusInfinite)

# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

add_subdirectory(Tgov1)
add_subdirectory(GASTPTI)
add_subdirectory(HYGOV)

install(FILES ResponseMode.hpp DESTINATION include/GridKit/Model/PhasorDynamics/Governor)

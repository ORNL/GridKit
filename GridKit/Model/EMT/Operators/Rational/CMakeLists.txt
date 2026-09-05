# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#]]

add_subdirectory(VectorFit)
add_subdirectory(StateSpace)

install(FILES Rational.hpp RationalMatrix.hpp RationalDataJSONParser.hpp
        DESTINATION include/GridKit/Model/EMT/Operators/Rational)

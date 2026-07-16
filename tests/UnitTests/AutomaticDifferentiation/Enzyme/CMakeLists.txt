# [[
#  Author(s):
#    - Nicholson Koukpaizan <koukpaizannk@ornl.gov>
#]]

add_executable(test_enzyme runEnzymeTests.cpp)
target_compile_options(test_enzyme PUBLIC -mllvm -enzyme-auto-sparsity=1)
target_link_libraries(test_enzyme GridKit::testing ClangEnzymeFlags)

add_test(NAME EnzymeTests COMMAND $<TARGET_FILE:test_enzyme>)

install(TARGETS test_enzyme 
        RUNTIME DESTINATION bin)

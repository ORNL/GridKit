# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(Dynamic)

if(TARGET solvers_opt_ipopt AND TARGET test_optimization_model)
  add_subdirectory(Optimization)
endif()

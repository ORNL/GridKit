gridkit_add_library(
  dense_vector
  SOURCES Vector.cpp
  HEADERS Vector.hpp
  LINK_LIBRARIES GridKit::utilities_logger GridKit::cpu_backend)

# # If no GPU is enabled link to dummy device backend
# if(NOT GridKit_ENABLE_GPU)
#   target_link_libraries(dense_vector PUBLIC GridKit::cpu_backend)
# endif(NOT GridKit_ENABLE_GPU)

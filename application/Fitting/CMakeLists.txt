# Fitting applications: rational approximation of sampled responses.
#
# Not yet wired into application/CMakeLists.txt; enable together with the
# GridKit::vector_fitting implementation.

add_executable(VectorFitting VectorFitting.cpp)
target_link_libraries(VectorFitting PUBLIC GridKit::vector_fitting
                                           GridKit::Utilities)
target_include_directories(
  VectorFitting PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

install(
  TARGETS VectorFitting
  EXPORT gridkit-targets
  RUNTIME)

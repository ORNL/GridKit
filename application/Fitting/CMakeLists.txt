# Fitting applications: rational approximation of sampled responses.
# fitting_model_json carries the shared rational-model JSON vocabulary.

add_library(fitting_model_json INTERFACE)
target_include_directories(
  fitting_model_json
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

add_executable(VectorFitting VectorFitting.cpp)
target_link_libraries(
  VectorFitting
  PRIVATE GridKit::vector_fitting GridKit::Utilities fitting_model_json)

install(
  TARGETS VectorFitting
  EXPORT gridkit-targets
  RUNTIME)

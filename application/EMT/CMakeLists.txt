add_executable(FrequencyResponse FrequencyResponse.cpp)
target_link_libraries(
  FrequencyResponse
  PUBLIC GridKit::emt_overhead GridKit::solvers_dyn GridKit::Utilities)
target_include_directories(
  FrequencyResponse
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

install(
  TARGETS FrequencyResponse
  EXPORT gridkit-targets
  RUNTIME)

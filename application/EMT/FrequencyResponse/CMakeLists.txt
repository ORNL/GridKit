# Frequency sweep engine shared by the EMT applications, and its driver.

add_library(emt_frequency_sweep INTERFACE)
target_include_directories(emt_frequency_sweep
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(emt_frequency_sweep INTERFACE GridKit::emt_overhead
                                                    GridKit::solvers_dyn)

add_executable(FrequencyResponse FrequencyResponse.cpp)
target_link_libraries(FrequencyResponse PRIVATE emt_frequency_sweep
                                                GridKit::Utilities)
target_include_directories(
  FrequencyResponse
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

install(
  TARGETS FrequencyResponse
  EXPORT gridkit-targets
  RUNTIME)

# UniversalLineModel: line description in, fitted line model out.
#
# Enabled from the parent directory once GridKit::vector_fitting is
# implemented.

add_executable(UniversalLineModel UniversalLineModel.cpp)
target_link_libraries(
  UniversalLineModel PRIVATE emt_frequency_sweep GridKit::vector_fitting
                             GridKit::Utilities)
target_include_directories(
  UniversalLineModel
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

install(
  TARGETS UniversalLineModel
  EXPORT gridkit-targets
  RUNTIME)

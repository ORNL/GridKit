# UniversalLineModel: line description in, fitted line model out.
# GridKit::Utilities carries the CliArgs interface; fitting_model_json the
# shared rational-model JSON vocabulary.

add_executable(UniversalLineModel UniversalLineModel.cpp)
target_link_libraries(
  UniversalLineModel
  PRIVATE emt_frequency_sweep
          GridKit::vector_fitting
          GridKit::Utilities
          fitting_model_json)

install(
  TARGETS UniversalLineModel
  EXPORT gridkit-targets
  RUNTIME)

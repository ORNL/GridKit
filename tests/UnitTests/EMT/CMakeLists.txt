if(GRIDKIT_ENABLE_ENZYME)
  add_executable(test_emt_models runModelTests.cpp)
  target_link_libraries(test_emt_models GridKit::emt_systemmodel GridKit::testing)
  target_include_directories(test_emt_models
                             PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
  add_test(NAME EMTModelTest COMMAND test_emt_models)
  install(TARGETS test_emt_models)
endif()

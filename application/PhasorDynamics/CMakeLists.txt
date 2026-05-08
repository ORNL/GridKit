add_executable(PDSim PDSim.cpp)
target_link_libraries(PDSim
    PUBLIC
    GridKit::phasor_dynamics_components
    GridKit::solvers_dyn
    GridKit::Utilities
    GridKit::testing)
target_include_directories(PDSim PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

find_package(OpenMP REQUIRED)
add_executable(ConCon ConCon.cpp CTGStudy.cpp)
target_link_libraries(ConCon
    PUBLIC
    OpenMP::OpenMP_CXX
    GridKit::phasor_dynamics_components
    GridKit::solvers_dyn
    GridKit::Utilities
    GridKit::testing)
target_include_directories(ConCon PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

install(TARGETS PDSim EXPORT gridkit-targets RUNTIME)
install(TARGETS ConCon EXPORT gridkit-targets RUNTIME)

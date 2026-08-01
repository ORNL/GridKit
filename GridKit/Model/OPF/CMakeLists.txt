set(opf_branch_headers Branch/Branch.hpp Branch/BranchData.hpp)
set(opf_bus_headers Bus/Bus.hpp Bus/BusData.hpp)
set(opf_generator_headers Generator/Generator.hpp Generator/GeneratorData.hpp)
set(opf_load_headers Load/Load.hpp Load/LoadData.hpp)
set(opf_shunt_headers Shunt/Shunt.hpp Shunt/ShuntData.hpp)

gridkit_add_library(
  model_opf
  SOURCES Branch/Branch.cpp
          Bus/Bus.cpp
          Generator/Generator.cpp
          Load/Load.cpp
          Shunt/Shunt.cpp
          System.cpp
          SystemData.cpp
  HEADERS Component.hpp
          ComponentConstraints.hpp
          ComponentVariables.hpp
          System.hpp
          SystemData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::model_state_data
    PUBLIC
    GridKit::dense_vector
    PUBLIC
    GridKit::sparse_matrix
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

target_sources(
  model_opf
  PRIVATE ${opf_branch_headers}
          ${opf_bus_headers}
          ${opf_generator_headers}
          ${opf_load_headers}
          ${opf_shunt_headers})

install(FILES ${opf_branch_headers} DESTINATION include/GridKit/Model/OPF/Branch)
install(FILES ${opf_bus_headers} DESTINATION include/GridKit/Model/OPF/Bus)
install(FILES ${opf_generator_headers} DESTINATION include/GridKit/Model/OPF/Generator)
install(FILES ${opf_load_headers} DESTINATION include/GridKit/Model/OPF/Load)
install(FILES ${opf_shunt_headers} DESTINATION include/GridKit/Model/OPF/Shunt)

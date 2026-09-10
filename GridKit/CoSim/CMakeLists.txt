gridkit_add_library(
  co_sim
  SOURCES CoSim.cpp
  HEADERS CoSim.hpp Message.hpp
  LINK_LIBRARIES
    PRIVATE
    cppzmq
  INCLUDE_DIRECTORIES
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

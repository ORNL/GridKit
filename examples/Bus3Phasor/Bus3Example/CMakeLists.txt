add_executable(bus3grid bus3phasor.cpp)
target_link_libraries(bus3grid 
	GRIDKIT::bus 
	SUNDIALS::sunlinsolklu
	SUNDIALS::core
	SUNDIALS::ida
	SUNDIALS::idas
	SUNDIALS::sunmatrixdense)
install(TARGETS bus3grid RUNTIME DESTINATION bin)
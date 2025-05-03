add_executable(gen2_example example.cpp)
target_link_libraries(gen2_example 
	GRIDKIT::bus 
	SUNDIALS::sunlinsolklu
	SUNDIALS::core
	SUNDIALS::ida
	SUNDIALS::idas
	SUNDIALS::sunmatrixdense)
install(TARGETS gen2_example RUNTIME DESTINATION bin)

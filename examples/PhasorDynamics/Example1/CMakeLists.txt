add_executable(phasordynamics_example1 example1.cpp)
target_link_libraries(phasordynamics_example1 
	GRIDKIT::bus 
	SUNDIALS::sunlinsolklu
	SUNDIALS::core
	SUNDIALS::ida
	SUNDIALS::idas
	SUNDIALS::sunmatrixdense)
install(TARGETS phasordynamics_example1 RUNTIME DESTINATION bin)

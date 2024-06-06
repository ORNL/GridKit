

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <math.h>

#include <ComponentLib/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>


/**
 * @brief Testing for the Distributed Generators outputs
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main(int argc, char const *argv[])
{
	
	ModelLib::DistributedGeneratorParameters<double, size_t> parms;
	//Parameters from MATLAB Microgrid code for first DG
	parms.wb_ = 2.0*M_PI*50.0;
	parms.wc_ = 31.41;
	parms.mp_ = 9.4e-5;
	parms.Vn_ = 380;
	parms.nq_ = 1.3e-3;
	parms.F_ = 0.75;
	parms.Kiv_ = 420.0;
	parms.Kpv_ = 0.1;
	parms.Kic_ = 20.0 * 1.0e3;
	parms.Kpc_ = 15.0;
	parms.Cf_ = 50.0e-6;
	parms.rLf_ = 0.1;
	parms.Lf_ = 1.35e-3;
	parms.rLc_ = 0.03;
	parms.Lc_ = 0.35e-3;

	ModelLib::DistributedGenerator<double, size_t> *dg = new ModelLib::DistributedGenerator<double, size_t>(0, parms, true);

	std::vector<double> t1(16,0.0);
	std::vector<double> t2{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};

	dg->allocate();

	dg->y() = t2;
	dg->yp() = t1;

	dg->evaluateResidual();

	std::cout << "Output: {";
	for (double i : dg->getResidual())
	{
		printf("%e ,", i);
	}
	std::cout << "}\n";

	//Generated from matlab code with same parameters and inputs
	std::vector<double> true_vec{3.141592277589793e+02,
     8.941907747838389e-01,
     1.846733023014284e+00,
     3.141592277589793e+02,
     1.014543000000000e+02,
    -1.507680000000000e+01,
     3.787993500000000e+02,
    -1.300000000000000e+00,
     2.899095146477517e+02,
     2.939138495559215e+02,
     1.507210571826699e+07,
     1.659799832843673e+07,
    -7.591593003913325e+03,
    -8.376991073310774e+03,
     3.337988298081817e+03,
     2.684419146397466e+03};

	std::cout << "Test the Relative Error\n";
	 for (size_t i = 0; i < true_vec.size(); i++)
	 {
		printf("%e ,\n", (true_vec[i] - dg->getResidual()[i]) / true_vec[i]);
	 }

	return 0;
}

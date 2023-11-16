

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <math.h>

#include <ComponentLib/PowerElectronicsComponents/DiscreteGenerator/DiscreteGenerator.hpp>


int main(int argc, char const *argv[])
{
	
	ModelLib::DiscreteGeneratorParameters<double, size_t> parms;
	//Parameters from MATLAB Microgrid code for first DG
	//refmp is need for reference input 
	parms.wb = 2.0*M_PI*50.0;
	parms.wc = 31.41;
	parms.mp = 9.4e-5;
	parms.refmp = 9.4e-5;
	parms.Vn = 380;
	parms.nq = 1.3e-3;
	parms.F = 0.75;
	parms.Kiv = 420.0;
	parms.Kpv = 0.1;
	parms.Kic = 20.0 * 1.0e3;
	parms.Kpc = 15.0;
	parms.Cf = 50.0e-6;
	parms.rLf = 0.1;
	parms.Lf = 1.35e-3;
	parms.rLc = 0.03;
	parms.Lc = 0.35e-3;

	ModelLib::DiscreteGenerator<double, size_t> *dg = new ModelLib::DiscreteGenerator<double, size_t>(0, parms);

	std::vector<double> t1(16,0.0);
	std::vector<double> t2{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};

	dg->allocate();

	dg->y() = t2;
	dg->yp() = t1;

	dg->evaluateResidual();

	std::cout << "Output: {";
	for (double i : dg->getResidual())
	{
		std::cout << i << ", ";
	}
	std::cout << "}\n";

	return 0;
}

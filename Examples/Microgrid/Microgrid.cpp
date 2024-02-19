

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>


#include <ComponentLib/PowerElectronicsComponents/DiscreteGenerator/DiscreteGenerator.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridLine/MicrogridLine.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridLoad/MicrogridLoad.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridBusDQ/MicrogridBusDQ.hpp>

#include <PowerElectronicsModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Solver/Dynamic/DynamicSolver.hpp>


int main(int argc, char const *argv[])
{
	double abstol = 1.0e-4;
	double reltol = 1.0e-4;
	bool usejac = true;

	//TODO:setup as named parameters
	//Create circuit model
	ModelLib::PowerElectronicsModel<double, size_t>* sysmodel = new ModelLib::PowerElectronicsModel<double, size_t>(reltol, abstol, usejac);

	//Modeled after the problem in the paper
	double RN = 1.0e4;

	//DG Params
	
	ModelLib::DiscreteGeneratorParameters<double, size_t> parms1;
	parms1.wb = 2.0*M_PI*50.0;
	parms1.wc = 31.41;
	parms1.mp = 9.4e-5;
	parms1.Vn = 380;
	parms1.nq = 1.3e-3;
	parms1.F = 0.75;
	parms1.Kiv = 420.0;
	parms1.Kpv = 0.1;
	parms1.Kic = 20.0 * 1.0e3;
	parms1.Kpc = 15.0;
	parms1.Cf = 50.0e-6;
	parms1.rLf = 0.1;
	parms1.Lf = 1.35e-3;
	parms1.rLc = 0.03;
	parms1.Lc = 0.35e-3;

	ModelLib::DiscreteGeneratorParameters<double, size_t> parms2;
	//Parameters from MATLAB Microgrid code for first DG
	parms2.wb = 2.0*M_PI*50.0;
	parms2.wc = 31.41;
	parms2.mp = 12.5e-5;
	parms2.Vn = 380;
	parms2.nq = 1.5e-3;
	parms2.F = 0.75;
	parms2.Kiv = 390.0;
	parms2.Kpv = 0.05;
	parms2.Kic = 16.0 * 1.0e3;
	parms2.Kpc = 10.5;
	parms2.Cf = 50.0e-6;
	parms2.rLf = 0.1;
	parms2.Lf = 1.35e-3;
	parms2.rLc = 0.03;
	parms2.Lc = 0.35e-3;

	//Line params
	double rline1 = 0.23;
	double Lline1 = 0.1 / (2.0 * M_PI * 50.0);

	double rline2 = 0.35;
	double Lline2 = 0.58 / (2.0 * M_PI * 50.0);

	double rline3 = 0.23;
	double Lline3 = 0.1 / (2.0 * M_PI * 50.0);

	//load parms
	double rload1 = 3.0;
	double Lload1 = 2.0 / (2.0 * M_PI * 50.0);

	double rload2 = 2.0;
	double Lload2 = 1.0 / (2.0 * M_PI * 50.0);


	//indexing sets
	size_t Nsize = 2;
	//							DGs	+			Lines +				Loads
	size_t vec_size_internals = 13*(2*Nsize) + (2 + 4*(Nsize - 1)) + 2*Nsize;
	//							\omegaref + BusDQ
	size_t vec_size_externals = 1 +	2*(2*Nsize);
	size_t dqbus1 = vec_size_externals + 1;
	size_t dqbus2 = vec_size_externals + 3;
	size_t dqbus3 = vec_size_externals + 5;
	size_t dqbus4 = vec_size_externals + 7;

	size_t vec_size_total = vec_size_internals + vec_size_externals;

	
	size_t indexv = 0;

	//dg 1
	ModelLib::DiscreteGenerator<double, size_t> *dg1 = new ModelLib::DiscreteGenerator<double, size_t>(0, parms1, true);
	//ref motor
	dg1->setExternalConnectionNodes(0,vec_size_internals);
	//outputs
	dg1->setExternalConnectionNodes(1,dqbus1);
	dg1->setExternalConnectionNodes(2,dqbus1 + 1);
	//internal connections
	for (size_t i = 0; i < 13; i++)
	{
		
		dg1->setExternalConnectionNodes(3 + i,indexv + i);
	}
	indexv += 13;
	sysmodel->addComponent(dg1);

	//dg 2
	ModelLib::DiscreteGenerator<double, size_t> *dg2 = new ModelLib::DiscreteGenerator<double, size_t>(1, parms1, false);
	//ref motor
	dg2->setExternalConnectionNodes(0,vec_size_internals);
	//outputs
	dg2->setExternalConnectionNodes(1,dqbus2);
	dg2->setExternalConnectionNodes(2,dqbus2 + 1);
	//internal connections
	for (size_t i = 0; i < 13; i++)
	{
		
		dg2->setExternalConnectionNodes(3 + i,indexv + i);
	}
	indexv += 13;
	sysmodel->addComponent(dg2);
	

	//dg 3
	ModelLib::DiscreteGenerator<double, size_t> *dg3 = new ModelLib::DiscreteGenerator<double, size_t>(2, parms2, false);
	//ref motor
	dg3->setExternalConnectionNodes(0,vec_size_internals);
	//outputs
	dg3->setExternalConnectionNodes(1,dqbus3);
	dg3->setExternalConnectionNodes(2,dqbus3 + 1);
	//internal connections
	for (size_t i = 0; i < 13; i++)
	{
		
		dg3->setExternalConnectionNodes(3 + i,indexv + i);
	}
	indexv += 13;
	sysmodel->addComponent(dg3);


	//dg 4
	ModelLib::DiscreteGenerator<double, size_t> *dg4 = new ModelLib::DiscreteGenerator<double, size_t>(3, parms2, false);
	//ref motor
	dg4->setExternalConnectionNodes(0,vec_size_internals);
	//outputs
	dg4->setExternalConnectionNodes(1,dqbus4);
	dg4->setExternalConnectionNodes(2,dqbus4 + 1);

	//internal connections
	for (size_t i = 0; i < 13; i++)
	{
		
		dg4->setExternalConnectionNodes(3 + i,indexv + i);
	}
	indexv += 13;
	sysmodel->addComponent(dg4);

	// Lines

	//line 1
	ModelLib::MicrogridLine<double, size_t> *l1 = new ModelLib::MicrogridLine<double, size_t>(4, rline1, Lline1);
	//ref motor
	l1->setExternalConnectionNodes(0,vec_size_internals);
	//input connections
	l1->setExternalConnectionNodes(1,dqbus1);
	l1->setExternalConnectionNodes(2,dqbus1 + 1);
	//output connections
	l1->setExternalConnectionNodes(3,dqbus2);
	l1->setExternalConnectionNodes(4,dqbus2 + 1);
	//internal connections
	for (size_t i = 0; i < 2; i++)
	{
		
		l1->setExternalConnectionNodes(5 + i,indexv + i);
	}
	indexv += 2;
	sysmodel->addComponent(l1);


	//line 2
	ModelLib::MicrogridLine<double, size_t> *l2 = new ModelLib::MicrogridLine<double, size_t>(5, rline2, Lline2);
	//ref motor
	l2->setExternalConnectionNodes(0,vec_size_internals);
	//input connections
	l2->setExternalConnectionNodes(1,dqbus2);
	l2->setExternalConnectionNodes(2,dqbus2 + 1);
	//output connections
	l2->setExternalConnectionNodes(3,dqbus3);
	l2->setExternalConnectionNodes(4,dqbus3 + 1);
	//internal connections
	for (size_t i = 0; i < 2; i++)
	{
		
		l2->setExternalConnectionNodes(5 + i,indexv + i);
	}
	indexv += 2;
	sysmodel->addComponent(l2);
	
	//line 3
	ModelLib::MicrogridLine<double, size_t> *l3 = new ModelLib::MicrogridLine<double, size_t>(6, rline3, Lline3);
	//ref motor
	l3->setExternalConnectionNodes(0,vec_size_internals);
	//input connections
	l3->setExternalConnectionNodes(1,dqbus3);
	l3->setExternalConnectionNodes(2,dqbus3 + 1);
	//output connections
	l3->setExternalConnectionNodes(3,dqbus4);
	l3->setExternalConnectionNodes(4,dqbus4 + 1);
	//internal connections
	for (size_t i = 0; i < 2; i++)
	{
		
		l3->setExternalConnectionNodes(5 + i,indexv + i);
	}
	indexv += 2;
	sysmodel->addComponent(l3);

	//  loads

	//load 1
	ModelLib::MicrogridLoad<double, size_t> *load1 = new ModelLib::MicrogridLoad<double, size_t>(7, rload1, Lload1);
	//ref motor
	load1->setExternalConnectionNodes(0,vec_size_internals);
	//input connections
	load1->setExternalConnectionNodes(1,dqbus1);
	load1->setExternalConnectionNodes(2,dqbus1 + 1);
	//internal connections
	for (size_t i = 0; i < 2; i++)
	{
		
		load1->setExternalConnectionNodes(3 + i,indexv + i);
	}
	indexv += 2;
	sysmodel->addComponent(load1);

	//load 2
	ModelLib::MicrogridLoad<double, size_t> *load2 = new ModelLib::MicrogridLoad<double, size_t>(8, rload2, Lload2);
	//ref motor
	load2->setExternalConnectionNodes(0,vec_size_internals);
	//input connections
	load2->setExternalConnectionNodes(1,dqbus3);
	load2->setExternalConnectionNodes(2,dqbus3 + 1);
	//internal connections
	for (size_t i = 0; i < 2; i++)
	{
		
		load2->setExternalConnectionNodes(3 + i,indexv + i);
	}
	indexv += 2;
	sysmodel->addComponent(load2);

	//Virtual PQ Buses
	ModelLib::MicrogridBusDQ<double, size_t> *bus1 = new ModelLib::MicrogridBusDQ<double, size_t>(9, RN);

	bus1->setExternalConnectionNodes(0,dqbus1);
	bus1->setExternalConnectionNodes(1,dqbus1 + 1);
	sysmodel->addComponent(bus1);

	ModelLib::MicrogridBusDQ<double, size_t> *bus2 = new ModelLib::MicrogridBusDQ<double, size_t>(10, RN);

	bus2->setExternalConnectionNodes(0,dqbus2);
	bus2->setExternalConnectionNodes(1,dqbus2 + 1);
	sysmodel->addComponent(bus2);

	ModelLib::MicrogridBusDQ<double, size_t> *bus3 = new ModelLib::MicrogridBusDQ<double, size_t>(11, RN);

	bus3->setExternalConnectionNodes(0,dqbus3);
	bus3->setExternalConnectionNodes(1,dqbus3 + 1);
	sysmodel->addComponent(bus3);

	ModelLib::MicrogridBusDQ<double, size_t> *bus4 = new ModelLib::MicrogridBusDQ<double, size_t>(12, RN);

	bus4->setExternalConnectionNodes(0,dqbus4);
	bus4->setExternalConnectionNodes(1,dqbus4 + 1);
	sysmodel->addComponent(bus4);

	sysmodel->allocate(vec_size_total);

	std::cout << sysmodel->y().size() << std::endl;
	std::cout << vec_size_internals << ", " << vec_size_externals << "\n";

	//Create Intial points for states
	for (size_t i = 0; i < vec_size_total; i++)
	{
		sysmodel->y()[i] = 0.0;
		sysmodel->yp()[i] = 0.0;
	}

	// Create Intial derivatives specifics generated in MATLAB
	//DGs
	for (size_t i = 0; i < 2; i++)
	{
		sysmodel->yp()[13*i + 3] = parms1.Vn;
		sysmodel->yp()[13*i + 5] = parms1.Kpv * parms1.Vn;
		sysmodel->yp()[13*i + 7] = (parms1.Kpc * parms1.Kpv * parms1.Vn) / parms1.Lf;
	}
	for (size_t i = 2; i < 4; i++)
	{
		sysmodel->yp()[13*i + 3] = parms2.Vn;
		sysmodel->yp()[13*i + 5] = parms2.Kpv * parms2.Vn;
		sysmodel->yp()[13*i + 7] = (parms2.Kpc * parms2.Kpv * parms2.Vn) / parms2.Lf;
	}

	//since the intial P_com = 0
	sysmodel->y()[62] = parms1.wb;

	

	sysmodel->initialize();
	sysmodel->evaluateResidual();

	std::vector<double>& fres = sysmodel->getResidual(); 
	std::cout << "Verify Intial Resisdual is Zero: {\n";
	for (size_t i = 0; i < fres.size(); i++)
	{
		printf("%u : %e \n", i, fres[i]);
	}
	std::cout << "}\n";

	sysmodel->updateTime(0.0, 1.0e-8);
	sysmodel->evaluateJacobian();
	std::cout << "Intial Jacobian with alpha:\n";
	// std::cout << sysmodel->getJacobian().frobnorm() << "\n";


	// //Create numerical integrator and configure it for the generator model
	// AnalysisManager::Sundials::Ida<double, size_t>* idas = new AnalysisManager::Sundials::Ida<double, size_t>(sysmodel);

    // double t_init  = 0.0;
    // double t_final = 1.0;

    // // setup simulation
    // idas->configureSimulation();
    // idas->getDefaultInitialCondition();
    // idas->initializeSimulation(t_init);

    // idas->runSimulation(t_final);

	// std::vector<double>& yfinial = sysmodel->y(); 

	// std::cout << "Final Vector y\n";
	// for (size_t i = 0; i < yfinial.size(); i++)
	// {
	// 	std::cout << yfinial[i] << "\n";
	// }

	return 0;
}

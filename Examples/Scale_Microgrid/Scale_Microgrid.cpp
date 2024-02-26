

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
	double abstol = 1.0e-8;
	double reltol = 1.0e-8;
	size_t max_step_amount = 3000;
	bool usejac = true;

	//TODO:setup as named parameters
	//Create circuit model
	ModelLib::PowerElectronicsModel<double, size_t>* sysmodel = new ModelLib::PowerElectronicsModel<double, size_t>(reltol, abstol, usejac, max_step_amount);

	
	//The amount of DG line load cobinations to generate for scale
	size_t Nsize = 2;

	//Modeled after the problem in the paper
	// Every Bus has the same virtual resistance. This is due to the numerical stability as mentioned in the paper.
	double RN = 1.0e4;

	//DG Params Vector
	//All DGs have the same set of parameters except for the first two.
	ModelLib::DiscreteGeneratorParameters<double, size_t> DG_parms1;
	DG_parms1.wb = 2.0*M_PI*50.0;
	DG_parms1.wc = 31.41;
	DG_parms1.mp = 9.4e-5;
	DG_parms1.Vn = 380.0;
	DG_parms1.nq = 1.3e-3;
	DG_parms1.F = 0.75;
	DG_parms1.Kiv = 420.0;
	DG_parms1.Kpv = 0.1;
	DG_parms1.Kic = 2.0e4;
	DG_parms1.Kpc = 15.0;
	DG_parms1.Cf = 5.0e-5;
	DG_parms1.rLf = 0.1;
	DG_parms1.Lf = 1.35e-3;
	DG_parms1.rLc = 0.03;
	DG_parms1.Lc = 0.35e-3;

	ModelLib::DiscreteGeneratorParameters<double, size_t> DG_parms2;
	DG_parms2.wb = 2.0*M_PI*50.0;
	DG_parms2.wc = 31.41;
	DG_parms2.mp = 12.5e-5;
	DG_parms2.Vn = 380.0;
	DG_parms2.nq = 1.5e-3;
	DG_parms2.F = 0.75;
	DG_parms2.Kiv = 390.0;
	DG_parms2.Kpv = 0.05;
	DG_parms2.Kic = 16.0e3;
	DG_parms2.Kpc = 10.5;
	DG_parms2.Cf = 50.0e-6;
	DG_parms2.rLf = 0.1;
	DG_parms2.Lf = 1.35e-3;
	DG_parms2.rLc = 0.03;
	DG_parms2.Lc = 0.35e-3;

	std::vector<ModelLib::DiscreteGeneratorParameters<double, size_t>> DGParams_list(2*Nsize, DG_parms2);

	DGParams_list[0] = DG_parms1;
	DGParams_list[1] = DG_parms1;

	//line vector params
	//Every odd line has the same parameters and every even line has the same parameters
	double rline1 = 0.23;
	double Lline1 = 0.1 / (2.0 * M_PI * 50.0);
	double rline2 = 0.35;
	double Lline2 = 0.58 / (2.0 * M_PI * 50.0);
	std::vector<double> rline_list(2*Nsize-1, 0.0);
	std::vector<double> Lline_list(2*Nsize-1, 0.0);
	for (size_t i = 0; i < rline_list.size(); i++)
	{
		rline_list[i] = (i % 2) ? rline1 : rline2;
		Lline_list[i] = (i % 2) ? Lline1 : Lline2;
	}
	

	//load parms
	//Only the first load has the same paramaters.
	double rload1 = 3.0;
	double Lload1 = 2.0 / (2.0 * M_PI * 50.0);
	double rload2 = 2.0;
	double Lload2 = 1.0 / (2.0 * M_PI * 50.0);

	std::vector<double> rload_list(Nsize, rload2);
	std::vector<double> Lload_list(Nsize, Lload2);
	rload_list[0] = rload1;
	Lload_list[0] = Lload1;

	//							DGs	+		- refframe	   Lines +				Loads
	size_t vec_size_internals = 13*(2*Nsize) - 1 + (2 + 4*(Nsize - 1)) + 2*Nsize;
	//							\omegaref + BusDQ
	size_t vec_size_externals = 1 +	2*(2*Nsize);

	std::vector<size_t> vdqbus_index(2*Nsize,0);
	vdqbus_index[0] = vec_size_internals + 1;
	for (size_t i = 1; i < vdqbus_index.size(); i++)
	{
		vdqbus_index[i] = vdqbus_index[i-1] + 2;
	}

	//Total size of the vector setup
	size_t vec_size_total = vec_size_internals + vec_size_externals;


	//Create the reference DG
	ModelLib::DiscreteGenerator<double, size_t> *dg_ref = new ModelLib::DiscreteGenerator<double, size_t>(0, DGParams_list[0], true);
	//ref motor
	dg_ref->setExternalConnectionNodes(0,vec_size_internals);
	//outputs
	dg_ref->setExternalConnectionNodes(1,vdqbus_index[0]);
	dg_ref->setExternalConnectionNodes(2,vdqbus_index[0] + 1);
	//"grounding" of the difference
	dg_ref->setExternalConnectionNodes(3,-1);
	//internal connections
	for (size_t i = 0; i < 12; i++)
	{
		
		dg_ref->setExternalConnectionNodes(4 + i,i);
	}
	sysmodel->addComponent(dg_ref);

	//Keep track of models and index location
	size_t indexv = 12;
	size_t model_id = 1;
	//Add all other DGs
	for (size_t i = 1; i < 2*Nsize; i++)
	{
		
		//current DG to add
		ModelLib::DiscreteGenerator<double, size_t> *dg = new ModelLib::DiscreteGenerator<double, size_t>(model_id, DGParams_list[i], false);
		//ref motor
		dg->setExternalConnectionNodes(0,vec_size_internals);
		//outputs
		dg->setExternalConnectionNodes(1,vdqbus_index[i]);
		dg->setExternalConnectionNodes(2,vdqbus_index[i] + 1);
		//internal connections
		for (size_t j = 0; j < 13; j++)
		{
			
			dg->setExternalConnectionNodes(3 + j,indexv + j);
		}
		indexv += 13;
		sysmodel->addComponent(dg);
		model_id++;
	}

	// Load all the Line compoenents
	for (size_t i = 0; i < 2*Nsize - 1; i++)
	{
		//line
		ModelLib::MicrogridLine<double, size_t> *line_model = new ModelLib::MicrogridLine<double, size_t>(model_id, rline_list[i], Lline_list[i]);
		//ref motor
		line_model->setExternalConnectionNodes(0,vec_size_internals);
		//input connections
		line_model->setExternalConnectionNodes(1,vdqbus_index[i]);
		line_model->setExternalConnectionNodes(2,vdqbus_index[i] + 1);
		//output connections
		line_model->setExternalConnectionNodes(3,vdqbus_index[i+1]);
		line_model->setExternalConnectionNodes(4,vdqbus_index[i+1] + 1);
		//internal connections
		for (size_t j = 0; j < 2; j++)
		{
			line_model->setExternalConnectionNodes(5 + j,indexv + j);
		}
		indexv += 2;
		sysmodel->addComponent(line_model);
		model_id++;
	}

	//  Load all the Load components
	for (size_t i = 0; i < Nsize; i++)
	{
		ModelLib::MicrogridLoad<double, size_t> *load_model = new ModelLib::MicrogridLoad<double, size_t>(model_id, rload_list[i], Lload_list[i]);
		//ref motor
		load_model->setExternalConnectionNodes(0,vec_size_internals);
		//input connections
		load_model->setExternalConnectionNodes(1,vdqbus_index[2*i]);
		load_model->setExternalConnectionNodes(2,vdqbus_index[2*i] + 1);
		//internal connections
		for (size_t j = 0; j < 2; j++)
		{
			
			load_model->setExternalConnectionNodes(3 + j,indexv + j);
		}
		indexv += 2;
		sysmodel->addComponent(load_model);
		model_id++;
	}

	//Add all the microgrid Virtual DQ Buses
	for (size_t i = 0; i < 2*Nsize; i++)
	{
		ModelLib::MicrogridBusDQ<double, size_t> *virDQbus_model = new ModelLib::MicrogridBusDQ<double, size_t>(
		model_id, RN);

		virDQbus_model->setExternalConnectionNodes(0, vdqbus_index[i]);
		virDQbus_model->setExternalConnectionNodes(1, vdqbus_index[i] + 1);
		sysmodel->addComponent(virDQbus_model);

		model_id++;
	}
	
	//allocate all the intial conditions
	sysmodel->allocate(vec_size_total);

	std::cout << sysmodel->y().size() << std::endl;
	std::cout << vec_size_internals << ", " << vec_size_externals << "\n";

	//Create Intial points for states. Every state is to specified to the zero intially
	for (size_t i = 0; i < vec_size_total; i++)
	{
		sysmodel->y()[i] = 0.0;
		sysmodel->yp()[i] = 0.0;
	}

	// Create Intial derivatives specifics generated in MATLAB
	for (size_t i = 0; i < 2*Nsize; i++)
	{
		sysmodel->yp()[13*i - 1 + 3] = DGParams_list[i].Vn;
		sysmodel->yp()[13*i - 1 + 5] = DGParams_list[i].Kpv * DGParams_list[i].Vn;
		sysmodel->yp()[13*i - 1 + 7] = (DGParams_list[i].Kpc * DGParams_list[i].Kpv * DGParams_list[i].Vn) / DGParams_list[i].Lf;
	}

	//since the intial P_com = 0, the set the intial vector to the reference frame
	sysmodel->y()[vec_size_internals] = DG_parms1.wb;

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


	//Create numerical integrator and configure it for the generator model
	AnalysisManager::Sundials::Ida<double, size_t>* idas = new AnalysisManager::Sundials::Ida<double, size_t>(sysmodel);

    double t_init  = 0.0;
    double t_final = 1.0;

    // setup simulation
    idas->configureSimulation();
    idas->getDefaultInitialCondition();
    idas->initializeSimulation(t_init);

    idas->runSimulation(t_final);

	std::vector<double>& yfinial = sysmodel->y(); 

	std::cout << "Final Vector y\n";
	for (size_t i = 0; i < yfinial.size(); i++)
	{
		std::cout << yfinial[i] << "\n";
	}

	if(Nsize == 2)
	{
		//Generate from MATLAB code ODE form with tolerances of 1e-14
		std::vector<double>true_vec{
		2.297543153595780e+04,
		1.275311524125022e+04,
		3.763060183116022e-02,
		-2.098153459325261e-02,
		1.848285659119097e-02,
		-1.563291404944864e-04,
		6.321941907011718e+01,
		-2.942264300846256e+01,
		3.634209302905854e+02,
		-2.668928293656362e-06,
		6.321941919221522e+01,
		-3.509200178595996e+01,
		-7.555954467454730e-03,
		2.297580486511343e+04,
		8.742028429066131e+03,
		3.710079564796484e-02,
		-1.421122598056797e-02,
		1.874079517807597e-02,
		-9.891304812687215e-05,
		6.232933298360234e+01,
		-1.796494061423331e+01,
		3.686353885026506e+02,
		3.465673854181523e-05,
		6.232933406188410e+01,
		-2.371564475187742e+01,
		-8.273939686941580e-02,
		1.727775042678524e+04,
		1.649365247247288e+04,
		3.116555157570849e-02,
		-2.985990066758010e-02,
		2.250012115906506e-02,
		-2.643873146501096e-04,
		4.861823510250247e+01,
		-4.088592755441309e+01,
		3.552597163751238e+02,
		-1.496407194199739e-04,
		4.861823504694532e+01,
		-4.642797132602495e+01,
		-8.445727984408551e-02,
		1.727723725566433e+04,
		9.182386962936238e+03,
		3.024959333190777e-02,
		-1.617250828202081e-02,
		2.318056864131751e-02,
		-1.295918667730514e-04,
		4.718938244522050e+01,
		-1.935782085675469e+01,
		3.662262287803608e+02,
		1.076423957830039e-04,
		4.718938116520511e+01,
		-2.507094256286497e+01,
		-1.881248349415025e+01,
		2.114714832305742e+01,
		4.329946674909793e+01,
		-3.037887936225145e+00,
		-4.487023117352992e+01,
		2.895883729832657e+01,
		8.199613345691378e+01,
		-5.623856502948122e+01,
		1.327498499660322e+02,
		-8.228065162347022e+01,
		3.119995747945993e+02,
		3.576922945168803e+02,
		-5.850795361581618e+00,
		3.641193316268954e+02,
		-8.846325267612976e+00,
		3.472146752739036e+02,
		-3.272400970143252e+01,
		3.604108939430972e+02,
		-3.492842627398574e+01
		};

		std::cout << "Test the Relative Error\n";
		for (size_t i = 0; i < true_vec.size(); i++)
		{
			printf("%u : %e ,\n", i, abs(true_vec[i] - yfinial[i]) / abs(true_vec[i]));
		}
	}

	return 0;
}

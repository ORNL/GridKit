


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>


#include <ComponentLib/PowerElectronicsComponents/DistributedGenerator/DistributedGenerator.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridLine/MicrogridLine.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridLoad/MicrogridLoad.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridBusDQ/MicrogridBusDQ.hpp>

#include <PowerElectronicsModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Solver/Dynamic/DynamicSolver.hpp>


int main(int argc, char const *argv[])
{
	///@todo Needs to be modified. Some components are small relative to others thus there error is high (or could be matlab vector issue)
    double abstol = 1.0e-8;
    double reltol = 1.0e-8;
    size_t max_step_amount = 3000;
    bool usejac = true;

    //Create model
    ModelLib::PowerElectronicsModel<double, size_t>* sysmodel = new ModelLib::PowerElectronicsModel<double, size_t>(reltol, abstol, usejac, max_step_amount);

    //Modeled after the problem in the paper
    double RN = 1.0e4;

    //DG Params
    ModelLib::DistributedGeneratorParameters<double, size_t> parms1;
    parms1.wb_ = 2.0*M_PI*50.0;
    parms1.wc_ = 31.41;
    parms1.mp_ = 9.4e-5;
    parms1.Vn_ = 380.0;
    parms1.nq_ = 1.3e-3;
    parms1.F_ = 0.75;
    parms1.Kiv_ = 420.0;
    parms1.Kpv_ = 0.1;
    parms1.Kic_ = 2.0e4;
    parms1.Kpc_ = 15.0;
    parms1.Cf_ = 5.0e-5;
    parms1.rLf_ = 0.1;
    parms1.Lf_ = 1.35e-3;
    parms1.rLc_ = 0.03;
    parms1.Lc_ = 0.35e-3;

    ModelLib::DistributedGeneratorParameters<double, size_t> parms2;
    //Parameters from MATLAB Microgrid code for first DG
    parms2.wb_ = 2.0*M_PI*50.0;
    parms2.wc_ = 31.41;
    parms2.mp_ = 12.5e-5;
    parms2.Vn_ = 380.0;
    parms2.nq_ = 1.5e-3;
    parms2.F_ = 0.75;
    parms2.Kiv_ = 390.0;
    parms2.Kpv_ = 0.05;
    parms2.Kic_ = 16.0e3;
    parms2.Kpc_ = 10.5;
    parms2.Cf_ = 50.0e-6;
    parms2.rLf_ = 0.1;
    parms2.Lf_ = 1.35e-3;
    parms2.rLc_ = 0.03;
    parms2.Lc_ = 0.35e-3;

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
    //							DGs	+		- refframe	   Lines +				Loads
    size_t vec_size_internals = 13*(2*Nsize) - 1 + (2 + 4*(Nsize - 1)) + 2*Nsize;
    //							\omegaref + BusDQ
    size_t vec_size_externals = 1 +	2*(2*Nsize);
    size_t dqbus1 = vec_size_internals + 1;
    size_t dqbus2 = vec_size_internals + 3;
    size_t dqbus3 = vec_size_internals + 5;
    size_t dqbus4 = vec_size_internals + 7;

    size_t vec_size_total = vec_size_internals + vec_size_externals;

    
    size_t indexv = 0;

    //dg 1
    ModelLib::DistributedGenerator<double, size_t> *dg1 = new ModelLib::DistributedGenerator<double, size_t>(0, parms1, true);
    //ref motor
    dg1->setExternalConnectionNodes(0,vec_size_internals);
    //outputs
    dg1->setExternalConnectionNodes(1,dqbus1);
    dg1->setExternalConnectionNodes(2,dqbus1 + 1);
    //"grounding" of the difference
    dg1->setExternalConnectionNodes(3,-1);
    //internal connections
    for (size_t i = 0; i < 12; i++)
    {
        
        dg1->setExternalConnectionNodes(4 + i,indexv + i);
    }
    indexv += 12;
    sysmodel->addComponent(dg1);

    //dg 2
    ModelLib::DistributedGenerator<double, size_t> *dg2 = new ModelLib::DistributedGenerator<double, size_t>(1, parms1, false);
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
    ModelLib::DistributedGenerator<double, size_t> *dg3 = new ModelLib::DistributedGenerator<double, size_t>(2, parms2, false);
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
    ModelLib::DistributedGenerator<double, size_t> *dg4 = new ModelLib::DistributedGenerator<double, size_t>(3, parms2, false);
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
    //DGs 1
    sysmodel->yp()[2] = parms1.Vn_;
    sysmodel->yp()[4] = parms1.Kpv_ * parms1.Vn_;
    sysmodel->yp()[6] = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
    sysmodel->yp()[12 + 3] = parms1.Vn_;
    sysmodel->yp()[12 + 5] = parms1.Kpv_ * parms1.Vn_;
    sysmodel->yp()[12 + 7] = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
    for (size_t i = 2; i < 4; i++)
    {
        sysmodel->yp()[13*i - 1 + 3] = parms2.Vn_;
        sysmodel->yp()[13*i - 1 + 5] = parms2.Kpv_ * parms2.Vn_;
        sysmodel->yp()[13*i - 1 + 7] = (parms2.Kpc_ * parms2.Kpv_ * parms2.Vn_) / parms2.Lf_;
    }

    //since the intial P_com = 0
    sysmodel->y()[vec_size_internals] = parms1.wb_;

    

    sysmodel->initialize();
    sysmodel->evaluateResidual();

    std::vector<double>& fres = sysmodel->getResidual(); 
    std::cout << "Verify Intial Resisdual is Zero: {\n";
    for (size_t i = 0; i < fres.size(); i++)
    {
        printf("%lu : %e \n", i, fres[i]);
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

    //Generate from MATLAB code ODE form with tolerances of 1e-12
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
    printf("%lu : %e ,\n", i, abs(true_vec[i] - yfinial[i]) / abs(true_vec[i]));
    }

	delete idas;
	delete sysmodel;

    return 0;
}

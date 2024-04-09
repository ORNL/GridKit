

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>

#include <ComponentLib/PowerElectronicsComponents/Capacitor/Capacitor.hpp>
#include <ComponentLib/PowerElectronicsComponents/Inductor/Inductor.hpp>
#include <ComponentLib/PowerElectronicsComponents/Resistor/Resistor.hpp>
#include <ComponentLib/PowerElectronicsComponents/VoltageSource/VoltageSource.hpp>

#include <PowerElectronicsModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Solver/Dynamic/DynamicSolver.hpp>


int main(int argc, char const *argv[])
{
    double abstol = 1.0e-8;
    double reltol = 1.0e-8;
    bool usejac = true;

    //TODO:setup as named parameters
    //Create circuit model
    ModelLib::PowerElectronicsModel<double, size_t>* sysmodel = new ModelLib::PowerElectronicsModel<double, size_t>(reltol, abstol, usejac);

    size_t idoff = 0;

    //RL circuit parameters
    double rinit = 1.0;
    double linit = 1.0;
    double vinit = 1.0;


    //inductor
    ModelLib::Inductor<double, size_t>* induct = new ModelLib::Inductor(idoff,linit);
    //Form index to node uid realations
    // input
    induct->setExternalConnectionNodes(0,1);
    //output
    induct->setExternalConnectionNodes(1,-1);
    //internal
    induct->setExternalConnectionNodes(2,2);
    //add component
    sysmodel->addComponent(induct);


    //resistor
    idoff++;
    ModelLib::Resistor<double, size_t>* resis = new ModelLib::Resistor(idoff, rinit);
    //Form index to node uid realations
    //input
    resis->setExternalConnectionNodes(0,0);
    //output
    resis->setExternalConnectionNodes(1,1);
    //add
    sysmodel->addComponent(resis);

    //voltage source
    idoff++;
    ModelLib::VoltageSource<double, size_t>* vsource = new ModelLib::VoltageSource(idoff, vinit);
    //Form index to node uid realations
    //input
    vsource->setExternalConnectionNodes(0,-1);
    //output
    vsource->setExternalConnectionNodes(1,0);
    //internal
    vsource->setExternalConnectionNodes(2,3);
    
    
    sysmodel->addComponent(vsource);

    sysmodel->allocate(4);

    std::cout << sysmodel->y().size() << std::endl;

    //Grounding for IDA. If no grounding then circuit is \mu > 1
    //v_0 (grounded)
    //Create Intial points
    sysmodel->y()[0] = vinit; //v_1
    sysmodel->y()[1] = vinit; // v_2
    sysmodel->y()[2] = 0.0; // i_L
    sysmodel->y()[3] = 0.0; // i_s

    sysmodel->yp()[0] = 0.0; // v'_1
    sysmodel->yp()[1] = 0.0; // v'_2
    sysmodel->yp()[2] = -vinit / linit; // i'_s
    sysmodel->yp()[3] = -vinit / linit; // i'_L
    
    
    sysmodel->initialize();
    sysmodel->evaluateResidual();

    std::cout << "Verify Intial Resisdual is Zero: {";
    for (double i : sysmodel->getResidual())
    {
        std::cout << i << ", ";
    }
    std::cout << "}\n";


    sysmodel->updateTime(0.0, 1.0);
    sysmodel->evaluateJacobian();
    std::cout << "Intial Jacobian with alpha = 1:\n";
    sysmodel->getJacobian().printMatrix();


    // Create numerical integrator and configure it for the generator model
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

    std::vector<double> yexact(4);

    //analytical solution to the circuit
    yexact[0] = vinit;
    yexact[2] = (vinit / rinit) * (exp(-(rinit / linit) * t_final) - 1.0);
    yexact[3] = yexact[2];
    yexact[1] = vinit + rinit * yexact[2];
    
    std::cout << "Element-wise Relative error at t=" << t_final << "\n";
    for (size_t i = 0; i < yfinial.size(); i++)
    {
        std::cout <<  abs((yfinial[i] - yexact[i]) / yexact[i]) << "\n";
    }

    return 0;
}

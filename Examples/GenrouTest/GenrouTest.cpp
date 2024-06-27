// Made by Paul Moon 6/7/2024


#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>

#include <ComponentLib/PowerFlow/Bus/BusSlack.hpp>
#include <ComponentLib/DynamicPhasor/SynchronousMachine/GENROUwS/GENROU.hpp>
#include <SystemModel.hpp>
#include <Solver/Dynamic/Ida.hpp>

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <Solver/Optimization/DynamicObjective.hpp>
#include <Solver/Optimization/DynamicConstraint.hpp>
#include <Utilities/Testing.hpp>


int main()
{
    using namespace ModelLib;
    using namespace AnalysisManager::Sundials;
    using namespace AnalysisManager;
    using namespace GridKit::Testing;
    

    // Create an infinite bus
    BaseBus<double, size_t>* bus = new BusSlack<double, size_t>(1.0946, 0.202); // roughly 1.072+0.22j

    // Attach a generator to that bus
    GENROU<double, size_t>* gen = new GENROU<double, size_t>(bus);

    // Create a system model
    SystemModel<double, size_t>* model = new SystemModel<double, size_t>();
    model->addBus(bus);
    model->addComponent(gen);

    // allocate model components
    model->allocate();

    std::cout << "Size: " << model->y().size() << std::endl;

    model->initialize();
    model->evaluateResidual();

    std::cout << "Verify Intial Residual is Zero: {";
    for (double i : model->getResidual())
    {
        std::cout << i << ", ";
    }
    std::cout << "}\n";


    //model->updateTime(0.0, 1.0);
    //model->evaluateJacobian();
    std::cout << "Intial Jacobian with alpha = 1:\n";
    //model->getJacobian().printMatrix();
    

    // Create numerical integrator and configure it for the generator model
    AnalysisManager::Sundials::Ida<double, size_t>* idas = new AnalysisManager::Sundials::Ida<double, size_t>(model);

    double t_init  = 0.0;
    double t_final = 2.0;
    double t_timestep = 0.0001;

    idas->configureSimulation();
    idas->getDefaultInitialCondition();
    idas->initializeSimulation(t_init);
    idas->runSimulation(t_final, 2000);

    delete idas;
    delete model;
    return 0;
}

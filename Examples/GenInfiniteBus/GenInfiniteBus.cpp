/*
 *
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov>.
 * LLNL-CODE-718378.
 * All rights reserved.
 *
 * This file is part of GridKit™. For details, see github.com/LLNL/GridKit
 * Please also read the LICENSE file.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the disclaimer below.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the disclaimer (as noted below) in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of the LLNS/LLNL nor the names of its contributors may
 *   be used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
 * SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISINGIN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Lawrence Livermore National Laboratory is operated by Lawrence Livermore
 * National Security, LLC, for the U.S. Department of Energy, National
 * Nuclear Security Administration under Contract DE-AC52-07NA27344.
 *
 * This document was prepared as an account of work sponsored by an agency
 * of the United States government. Neither the United States government nor
 * Lawrence Livermore National Security, LLC, nor any of their employees
 * makes any warranty, expressed or implied, or assumes any legal liability
 * or responsibility for the accuracy, completeness, or usefulness of any
 * information, apparatus, product, or process disclosed, or represents that
 * its use would not infringe privately owned rights. Reference herein to
 * any specific commercial product, process, or service by trade name,
 * trademark, manufacturer, or otherwise does not necessarily constitute or
 * imply its endorsement, recommendation, or favoring by the United States
 * government or Lawrence Livermore National Security, LLC. The views and
 * opinions of authors expressed herein do not necessarily state or reflect
 * those of the United States government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */


#include <iostream>
#include <iomanip>

#include <ComponentLib/Bus/BusSlack.hpp>
#include <ComponentLib/Generator4/Generator4.hpp>
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
    BaseBus<double, size_t>* bus = new BusSlack<double, size_t>(1.0, 0.0);

    // Attach a generator to that bus
    Generator4<double, size_t>* gen = new Generator4<double, size_t>(bus);

    // Create a system model
    SystemModel<double, size_t>* model = new SystemModel<double, size_t>();
    model->addBus(bus);
    model->addComponent(gen);

    // allocate model components
    model->allocate();
    

    // Create numerical integrator and configure it for the generator model
    Ida<double, size_t>* idas = new Ida<double, size_t>(model);

    double t_init  = 0.0;
    double t_final = 15.0;

    // setup simulation
    idas->configureSimulation();
    idas->configureAdjoint();
    idas->getDefaultInitialCondition();
    idas->initializeSimulation(t_init);
    idas->configureQuadrature();
    idas->initializeQuadrature();

    double t_fault = 0.1;
    double t_clear = 0.1;
    idas->runSimulation(t_fault);
    idas->saveInitialCondition();
    // create initial condition after a fault
    {
        idas->getSavedInitialCondition();
        idas->initializeSimulation(t_init);
        gen->V() = 0.0;
        idas->runSimulation(t_clear, 20);
        gen->V() = 1.0;
        idas->saveInitialCondition();
    }

    // Set integration time for dynamic constrained optimization
    idas->setIntegrationTime(t_init, t_final, 100);

    // Guess initial values of optimization parameters
    double Pm = 1.0;
    double Ef = 1.45;

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = IpoptApplicationFactory();

    // Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status;
    status = ipoptApp->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << "\n\n*** Initialization failed! ***\n\n";
        return (int) status;
    }

    // Set solver tolerance
    const double tol = 1e-4; 

    // Configure Ipopt application
    ipoptApp->Options()->SetStringValue("hessian_approximation", "limited-memory");
    ipoptApp->Options()->SetNumericValue("tol", tol);
    ipoptApp->Options()->SetIntegerValue("print_level", 0);

    // Create dynamic objective interface to Ipopt solver
    Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicObjectiveInterface =
        new IpoptInterface::DynamicObjective<double, size_t>(idas);

    // Initialize the problem
    model->param()[0] = Pm;
    model->param()[1] = Ef;

    // Solve the problem
    status = ipoptApp->OptimizeTNLP(ipoptDynamicObjectiveInterface);
    std::cout << "\n\nProblem formulated as dynamic objective optimiztion ...\n";

    if (status == Ipopt::Solve_Succeeded) {
        // Print result
        std::cout << "\nSucess:\n The problem solved in "
                  << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
                  << " Optimal value of Pm = " << model->param()[0] << "\n"
                  << " Optimal value of Ef = " << model->param()[1] << "\n"
                  << " The final value of the objective function G(Pm,Ef) = "
                  << ipoptApp->Statistics()->FinalObjective() << "\n\n";
    }

    // Create dynamic constraint interface to Ipopt solver
    Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicConstraintInterface =
        new IpoptInterface::DynamicConstraint<double, size_t>(idas);

    // Store dynamic objective optimization results
    double* results  = new double[model->size_opt()];
    for(unsigned i=0; i <model->size_opt(); ++i)
    {
        results[i] = model->param()[i];
    }

    // Initialize the problem
    model->param()[0] = Pm;
    model->param()[1] = Ef;

    // Solve the problem
    status = ipoptApp->OptimizeTNLP(ipoptDynamicConstraintInterface);
    std::cout << "\n\nProblem formulated as dynamic constraint optimiztion ...\n";

    if (status == Ipopt::Solve_Succeeded) {
        // Print result
        std::cout << "\nSucess:\n The problem solved in "
                  << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
                  << " Optimal value of Pm = " << model->param()[0] << "\n"
                  << " Optimal value of Ef = " << model->param()[1] << "\n"
                  << " The final value of the objective function G(Pm,Ef) = "
                  << ipoptApp->Statistics()->FinalObjective() << "\n\n";
    }

    // Compare results of the two optimization methods
    int retval = 0;
    for(unsigned i=0; i <model->size_opt(); ++i)
    {
        if(!isEqual(results[i], model->param()[i], 100*tol))
            --retval; 
    }

    if(retval < 0)
    {
        std::cout << "The two results differ beyond solver tolerance!\n";
    }

    delete [] results;
    delete idas;
    delete model;
    return 0;
}

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

#include <ComponentLib/PowerFlow/Bus/BusSlack.hpp>
#include <ComponentLib/PowerFlow/Generator4/Generator4.hpp>
#include <SystemModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Utilities/Testing.hpp>

/*
 * Compute gradient of an objective function expressed as an integral over
 * system trajectory. The gradient is computed numerically and using
 * adjoint sensitivity analysis.
 *
 * The test case is a 4th order generator connected to an infinite bus.
 * The objective function is total frequency deviation computed over
 * system trajectory after generator short circuit fault.
 *
 */
int main()
{
    using namespace ModelLib;
    using namespace AnalysisManager::Sundials;
    using namespace AnalysisManager;
    using namespace GridKit::Testing;

    // Create an infinite bus
    BaseBus<double, size_t>* bus = new BusSlack<double, size_t>(1.0, 0.0);

    // Attach a generator to that bus
    Generator4<double, size_t>* gen = new Generator4<double, size_t>(bus, 0.8, 0.3);

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


    idas->runSimulation(0.1, 2);
    idas->saveInitialCondition();

    // create initial condition after a fault
    {
        idas->getSavedInitialCondition();
        idas->initializeSimulation(t_init);
        gen->V() = 0.0;
        idas->runSimulation(0.1, 20);
        gen->V() = 1.0;
        idas->saveInitialCondition();
    }

    // Get pointer the objective function
    const double* Q = idas->getIntegral();

    // Compute the objective function as an integral over the system trajectory
    idas->getSavedInitialCondition();
    idas->initializeSimulation(t_init);
    idas->initializeQuadrature();
    idas->runSimulationQuadrature(t_final, 100);

    std::cout << "\n\nCost of computing objective function:\n\n";
    idas->printFinalStats();

    const double g1 = Q[0];
    const double eps = 2e-3;

    // Compute gradient of the objective function numerically
    std::vector<double> dGdp(model->size_opt());

    for (unsigned i=0; i<model->size_opt(); ++i)
    {
      model->param()[i] += eps;
      idas->getSavedInitialCondition();
      idas->initializeSimulation(t_init);
      idas->initializeQuadrature();
      idas->runSimulationQuadrature(t_final,100);

      std::cout << "\n\nCost of computing derivative with respect to parameter "
                << i << ":\n\n";
      idas->printFinalStats();
      double g2 = Q[0];

      // restore parameter to original value
      model->param()[i] -= eps;

      // Evaluate dG/dp numerically
      dGdp[i] = (g2 - g1)/eps;
    }

    // Compute gradient of the objective function using adjoint method
    idas->initializeAdjoint();
    idas->getSavedInitialCondition();
    idas->initializeSimulation(t_init);
    idas->initializeQuadrature();
    idas->runForwardSimulation(t_final, 100);

    std::cout << "\n\nCost of forward simulation for adjoint\n"
              << "sensitivity analysis:\n\n";
    idas->printFinalStats();

    idas->initializeBackwardSimulation(t_final);
    idas->runBackwardSimulation(t_init);

    std::cout << "\n\nCost of adjoint sensitivity analysis:\n\n";
    idas->printFinalStats();

    // Compare results
    int retval = 0;
    std::cout << "\n\nComparison of numerical and adjoint results:\n\n";
    double* neg_dGdp = idas->getAdjointIntegral();
    for (unsigned i=0; i<model->size_opt(); ++i)
    {
        std::cout << "dG/dp" << i << " (numerical) = " <<      dGdp[i] << "\n";
        std::cout << "dG/dp" << i << " (adjoint)   = " << -neg_dGdp[i] << "\n\n";
        if(!isEqual(dGdp[i], -neg_dGdp[i], 10*eps))
            --retval; 
    }

    if(retval < 0)
    {
        std::cout << "The two results differ beyond solver tolerance!\n";
    }

    return retval;
}

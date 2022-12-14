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

/**
 * @file Grid3BusSys.cpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 * 
 * Simple 3-bus grid example. Two models are tested here -- a hard-wired model
 * and a model assembled using GridKit's system composer.
 * 
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include <ComponentLib/MiniGrid/MiniGrid.hpp>
#include <ComponentLib/Bus/BusFactory.hpp>
#include <ComponentLib/Generator/GeneratorFactory.hpp>
#include <ComponentLib/Branch/Branch.hpp>
#include <ComponentLib/Load/Load.hpp>
#include <SystemSteadyStateModel.hpp>
#include <Solver/SteadyState/Kinsol.hpp>
#include <PowerSystemData.hpp>

#include <Utilities/Testing.hpp>


int main()
{
    using namespace ModelLib;
    using namespace AnalysisManager::Sundials;
    using namespace AnalysisManager;
    using namespace GridKit::Testing;
    using namespace GridKit::PowerSystemData;

    std::cout << "\nSolving power flow for a 3-bus monolithic model ...\n\n";
    // Create a 3-bus model
    MiniGrid<double, size_t>* model = new MiniGrid<double, size_t>();

    // allocate model
    model->allocate();
    std::cout << "Model size: " << model->size() << "\n\n";

    // Create numerical solver and attach the model to it.
    // Here we use Kinsol solver from SUNDIALS library
    Kinsol<double, size_t>* kinsol = new Kinsol<double, size_t>(model);

    // setup simulation
    kinsol->configureSimulation();
    // initialize simulation with default initial guess V=1, theta=0
    kinsol->getDefaultInitialCondition();
    // Compute solution
    kinsol->runSimulation();
    // Print solution
    double th2 = model->th2() * 180.0/M_PI; 
    double V2  = model->V2();
    double th3 = model->th3() * 180.0/M_PI; 
    std::cout << "Solution:\n";
    std::cout << "  theta2 = " << th2 << " deg,  expected = " << " -4.87979 deg\n";
    std::cout << "  V2     = " << V2  << " p.u., expected = " << "  1.08281 p.u.\n";
    std::cout << "  theta3 = " << th3 << " deg,  expected = " << "  1.46241 deg\n\n";

    // Print solver performance statistics
    kinsol->printFinalStats();

    int retval1 = 0;
    retval1 += isEqual(th2, -4.878, 1e-4);
    retval1 += isEqual(V2,   1.096, 1e-4);
    retval1 += isEqual(th3,  1.491, 1e-4);

    if(retval1 == 0)
        std::cout << "\nSucess!\n\n\n";

    // Delete solver and model
    delete kinsol; kinsol = nullptr;
    delete model;  model  = nullptr;

    std::cout << "Solving same problem, but assembled from components ...\n\n";

    // First, create an empty system model
    SystemSteadyStateModel<double, size_t>* sysmodel = new SystemSteadyStateModel<double, size_t>();

    // Next create and add buses ...
    // Create a slack bus, fix V=1, theta=0, bus ID = 1" << std::endl;
    BusData<double, size_t> bd1;
    bd1.bus_i = 1; bd1.type = 3; bd1.Vm = 1.0; bd1.Va = 0.0;
    auto* bus1 = BusFactory<double, size_t>::create(bd1);
    // BaseBus<double, size_t>* bus1 = new BusSlack<double, size_t>(bd1);
    sysmodel->addBus(bus1);
    
    // Create a PQ bus, initialize V=1, theta=0, bus ID = 2" << std::endl;
    BusData<double, size_t> bd2;
    bd2.bus_i = 2; bd2.type = 1; bd2.Vm = 1.0; bd2.Va = 0.0;
    auto* bus2 = BusFactory<double, size_t>::create(bd2);
    // BaseBus<double, size_t>* bus2 = new BusPQ<double, size_t>(bd2);
    sysmodel->addBus(bus2);

    std::cout << "// Create a PV bus, fix V=1.1, initialize theta=0, and set power injection Pg=2" << std::endl;
    BusData<double, size_t> bd3;
    bd3.bus_i = 3; bd3.type = 2; bd3.Vm = 1.1; bd3.Va = 0.0;
    auto* bus3 = BusFactory<double, size_t>::create(bd3);
    // BaseBus<double, size_t>* bus3 = new BusPV<double, size_t>(1.1, 0.0, 2.0);
    sysmodel->addBus(bus3);

    // Create and add generators ...
    std::cout << "// Set power injection to 2.0 for bus3" << std::endl;
    GenData<double, size_t> gd3;
    gd3.Pg = 2.0; gd3.bus = 3;
    auto* gen3 = GeneratorFactory<double, size_t>::create(bus3, gd3);
    sysmodel->addComponent(gen3);

    // Create and add branches ...
    Branch<double, size_t>* branch12 = new Branch<double, size_t>(bus1, bus2);
    branch12->setX(1.0/10.0);
    sysmodel->addComponent(branch12);
    Branch<double, size_t>* branch13 = new Branch<double, size_t>(bus1, bus3);
    branch13->setX(1.0/15.0);
    sysmodel->addComponent(branch13);
    Branch<double, size_t>* branch23 = new Branch<double, size_t>(bus2, bus3);
    branch23->setX(1.0/12.0);
    sysmodel->addComponent(branch23);

    // Create and add loads ...
    Load<double, size_t>* load1 = new Load<double, size_t>(bus1, 2.0, 0.0);
    sysmodel->addComponent(load1);
    Load<double, size_t>* load2 = new Load<double, size_t>(bus2, 2.5, -0.8);
    sysmodel->addComponent(load2);

    // allocate model
    sysmodel->allocate();
    std::cout << "Model size: " << sysmodel->size() << "\n\n";

    // Create numerical solver and attach the model to it.
    // Here we use Kinsol solver from SUNDIALS library
    kinsol = new Kinsol<double, size_t>(sysmodel);

    // setup simulation
    kinsol->configureSimulation();
    // initialize simulation with default initial guess 
    kinsol->getDefaultInitialCondition();
    // Compute solution
    kinsol->runSimulation();
    // Print solution
    th2 = bus2->theta() * 180.0/M_PI; 
    V2  = bus2->V();
    th3 = bus3->theta() * 180.0/M_PI; 
    std::cout << "Solution:\n";
    std::cout << "  theta2 = " << th2 << " deg,  expected = " << " -4.87979 deg\n";
    std::cout << "  V2     = " << V2  << " p.u., expected = " << "  1.08281 p.u.\n";
    std::cout << "  theta3 = " << th3 << " deg,  expected = " << "  1.46241 deg\n\n";

    // Print solver performance statistics
    kinsol->printFinalStats();

    int retval2 = 0;
    retval2 += isEqual(th2, -4.878, 1e-4);
    retval2 += isEqual(V2,   1.096, 1e-4);
    retval2 += isEqual(th3,  1.491, 1e-4);

    if(retval2 == 0)
        std::cout << "\nSucess!\n\n\n";


    // Delete solver and model
    delete kinsol;
    delete sysmodel;
    return retval1 + retval2;
}

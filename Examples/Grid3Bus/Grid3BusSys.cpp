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
 * @author Reid Gomillion <rjg18@vt.edu>
 * 
 * Simple 3-bus grid example. Two models are tested here -- a hard-wired model
 * and a model assembled using GridKit's system composer.
 * 
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>

#include <ComponentLib/MiniGrid/MiniGrid.hpp>
#include <ComponentLib/Bus/BusFactory.hpp>
#include <ComponentLib/Generator/GeneratorFactory.hpp>
#include <ComponentLib/Branch/Branch.hpp>
#include <ComponentLib/Load/Load.hpp>
#include <SystemSteadyStateModel.hpp>
#include <Solver/SteadyState/Kinsol.hpp>
#include <PowerSystemData.hpp>

#include <Utilities/FileIO.hpp>
#include <Utilities/Testing.hpp>

static const std::string BUS3_DATA_STRING = R"(
function mpc = case5
% Created by Reid Gomillion

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
    1	3	2.0	0.0	0	0	0	1.0	0.0	0	0	0	0.0;
    2	1	2.5	-0.8	0	0	0	1.0	0.0	0	0	0	0.0;
    3	2	0	0	0	0	0	1.1	0.0	0	0	0	0.0;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    3	2.0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    1	2	0	0.1	0	0	0	0	0	0	0	0	0;
    1	3	0	0.0666666	0	0	0	0	0	0	0	0	0;
    2	3	0	0.0833333	0	0	0	0	0	0	0	0	0;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
    2	0	0	3	0   14	0;
    2	0	0	3	0   15	0;
    2	0	0	3	0   30	0;
];

)";


using namespace ModelLib;
using namespace AnalysisManager::Sundials;
using namespace AnalysisManager;
using namespace GridKit::Testing;
using namespace GridKit::PowerSystemData;

constexpr double theta2_ref = -4.87979; // [deg]
constexpr double V2_ref     =  1.08281; // [p.u.]
constexpr double theta3_ref =  1.46241; // [deg]

/**
 * Testing the monlithic case via the class MiniGrid
 * @return returns 0 if pass o.w. fails
*/
int monolithicCase()
{
     std::cout << "\nSolving power flow for a 3-bus monolithic model ...\n\n";
    // Create a 3-bus model
    MiniGrid<double, size_t>* model = new MiniGrid<double, size_t>();

    // allocate model
    model->allocate();
    model->initialize();
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
    std::cout << "  theta2 = " << th2 << " deg,  expected = " << theta2_ref << " deg\n";
    std::cout << "  V2     = " << V2  << " p.u., expected = " << V2_ref     << " p.u.\n";
    std::cout << "  theta3 = " << th3 << " deg,  expected = " << theta3_ref << " deg\n\n";

    // Print solver performance statistics
    kinsol->printFinalStats();

    int retval1 = 0;
    retval1 += !isEqual(th2, theta2_ref, 1e-4);
    retval1 += !isEqual(V2,  V2_ref    , 1e-4);
    retval1 += !isEqual(th3, theta3_ref, 1e-4);

    if(retval1 == 0)
        std::cout << "\nSuccess!\n\n\n";
    else
        std::cout << "\nFailed!\n\n\n";

    // Delete solver and model
    delete kinsol; kinsol = nullptr;
    delete model;  model  = nullptr;
    return retval1;
}

/**
 * Run the Testing case for parser setup
 * @return returns 0 if pass o.w. fail
*/
int parserCase()
{
    std::cout << "Solving same problem, but assembled from components via a parser ...\n\n";

    //Data File Reading
    GridKit::PowerSystemData::SystemModelData<double, size_t> mp;

    std::istringstream iss(BUS3_DATA_STRING);
    GridKit::readMatPower(mp, iss);

    
    // Create an empty system model and pass the system model data to it
    SystemSteadyStateModel<double, size_t>* sysmodel = new SystemSteadyStateModel<double, size_t>(mp);

    // allocate model
    sysmodel->allocate();
    sysmodel->initialize();
    std::cout << "Model size: " << sysmodel->size() << "\n\n";

    // Create numerical solver and attach the model to it.
    // Here we use Kinsol solver from SUNDIALS library
    Kinsol<double, size_t>* kinsol = new Kinsol<double, size_t>(sysmodel);

    // setup simulation
    kinsol->configureSimulation();
    // initialize simulation with default initial guess 
    kinsol->getDefaultInitialCondition();
    // Compute solution
    kinsol->runSimulation();
    // Print solution
    double th2 = sysmodel->getBus(2)->theta() * 180.0/M_PI;
    double V2 =  sysmodel->getBus(2)->V();
    double th3 = sysmodel->getBus(3)->theta() * 180.0/M_PI;


    std::cout << "Solution:\n";
    std::cout << "  theta2 = " << th2 << " deg,  expected = " << theta2_ref << " deg\n";
    std::cout << "  V2     = " << V2  << " p.u., expected = " << V2_ref     << " p.u.\n";
    std::cout << "  theta3 = " << th3 << " deg,  expected = " << theta3_ref << " deg\n\n";

    // Print solver performance statistics
    kinsol->printFinalStats();

    int retval2 = 0;
    retval2 += !isEqual(th2, theta2_ref, 1e-4);
    retval2 += !isEqual(V2,  V2_ref    , 1e-4);
    retval2 += !isEqual(th3, theta3_ref, 1e-4);

    if(retval2 == 0)
        std::cout << "\nSuccess!\n\n\n";
    else
        std::cout << "\nFailed!\n\n\n";


    // Delete solver and model
    
    delete kinsol; kinsol = nullptr;
    delete sysmodel;  sysmodel  = nullptr;
    
    return retval2;
}

/**
 * Hardwired Test Case
 * @return 0 if pass otherwise fails
*/
int hardwiredCase()
{
    std::cout << "Solving same problem, but assembled from components manually ...\n\n";

    // First, create an empty system model
    SystemSteadyStateModel<double, size_t>* sysmodel = new SystemSteadyStateModel<double, size_t>();

    // Next create and add buses ...
    // Create a slack bus, fix V=1, theta=0, bus ID = 1
    BusData<double, size_t> bd1;
    bd1.bus_i = 1; bd1.type = 3; bd1.Vm = 1.0; bd1.Va = 0.0;
    auto* bus1 = BusFactory<double, size_t>::create(bd1);
    sysmodel->addBus(bus1);
    
    //Create a PQ bus, initialize V=1, theta=0, bus ID = 2
    BusData<double, size_t> bd2;
    bd2.bus_i = 2; bd2.type = 1; bd2.Vm = 1.0; bd2.Va = 0.0;
    auto* bus2 = BusFactory<double, size_t>::create(bd2);
    sysmodel->addBus(bus2);

    // Create a PV bus, fix V=1.1, initialize theta=0, and set power injection Pg=2
    BusData<double, size_t> bd3;
    bd3.bus_i = 3; bd3.type = 2; bd3.Vm = 1.1; bd3.Va = 0.0;
    auto* bus3 = BusFactory<double, size_t>::create(bd3);
    sysmodel->addBus(bus3);

    // Create and add generators ...
    // Create and add slack generator connected to bus1
    GenData<double, size_t> gd1;
    gd1.bus = 1;
    auto* gen1 = GeneratorFactory<double, size_t>::create(sysmodel->getBus(gd1.bus), gd1);
    sysmodel->addComponent(gen1);

    // Create and add PV generator connected to bus3
    GenData<double, size_t> gd3;
    gd3.Pg = 2.0; gd3.bus = 3;
    auto* gen3 = GeneratorFactory<double, size_t>::create(sysmodel->getBus(gd3.bus), gd3);
    sysmodel->addComponent(gen3);

    // Create and add branches ...
    // Branch 1-2
    BranchData<double, size_t> brd12;
    brd12.fbus = 1; brd12.tbus = 2; brd12.x = 1.0/10.0; brd12.r = 0.0; brd12.b = 0.0;
    Branch<double, size_t>* branch12 = new Branch<double, size_t>(sysmodel->getBus(brd12.fbus), sysmodel->getBus(brd12.tbus), brd12);
    sysmodel->addComponent(branch12);

    // Branch 1-3
    BranchData<double, size_t> brd13;
    brd13.fbus = 1; brd13.tbus = 3; brd13.x = 1.0/15.0; brd13.r = 0.0; brd13.b = 0.0;
    Branch<double, size_t>* branch13 = new Branch<double, size_t>(sysmodel->getBus(brd13.fbus), sysmodel->getBus(brd13.tbus), brd13);
    sysmodel->addComponent(branch13);

    // Branch 2-3
    BranchData<double, size_t> brd23;
    brd23.fbus = 2; brd23.tbus = 3; brd23.x = 1.0/12.0; brd23.r = 0.0; brd23.b = 0.0;
    Branch<double, size_t>* branch23 = new Branch<double, size_t>(sysmodel->getBus(brd23.fbus), sysmodel->getBus(brd23.tbus), brd23);
    sysmodel->addComponent(branch23);
    

    // Create and add loads ...
    // Load on bus1
    LoadData<double, size_t> ld1;
    ld1.bus_i = 1; ld1.Pd = 2.0; ld1.Qd = 0.0;
    Load<double, size_t>* load1 = new Load<double, size_t>(sysmodel->getBus(ld1.bus_i), ld1);
    sysmodel->addComponent(load1);

    // Load on bus2
    LoadData<double, size_t> ld2;
    ld2.bus_i = 2; ld2.Pd = 2.5; ld2.Qd = -0.8;
    Load<double, size_t>* load2 = new Load<double, size_t>(sysmodel->getBus(ld2.bus_i), ld2);
    sysmodel->addComponent(load2);

    // allocate model
    sysmodel->allocate();
    sysmodel->initialize();
    std::cout << "Model size: " << sysmodel->size() << "\n\n";

    // Create numerical solver and attach the model to it.
    // Here we use Kinsol solver from SUNDIALS library
    Kinsol<double, size_t>* kinsol = new Kinsol<double, size_t>(sysmodel);

    // setup simulation
    kinsol->configureSimulation();
    // initialize simulation with default initial guess 
    kinsol->getDefaultInitialCondition();
    // Compute solution
    kinsol->runSimulation();
    // Print solution
    double th2 = bus2->theta() * 180.0/M_PI; 
    double V2  = bus2->V();
    double th3 = bus3->theta() * 180.0/M_PI;


    std::cout << "Solution:\n";
    std::cout << "  theta2 = " << th2 << " deg,  expected = " << theta2_ref << " deg\n";
    std::cout << "  V2     = " << V2  << " p.u., expected = " << V2_ref     << " p.u.\n";
    std::cout << "  theta3 = " << th3 << " deg,  expected = " << theta3_ref << " deg\n\n";

    // Print solver performance statistics
    kinsol->printFinalStats();

    int retval2 = 0;
    retval2 += !isEqual(th2, theta2_ref, 1e-4);
    retval2 += !isEqual(V2,  V2_ref    , 1e-4);
    retval2 += !isEqual(th3, theta3_ref, 1e-4);

    if(retval2 == 0)
        std::cout << "\nSuccess!\n\n\n";
    else
        std::cout << "\nFailed!\n\n\n";


    // Delete solver and model
    delete kinsol; kinsol = nullptr;
    delete sysmodel;  sysmodel  = nullptr;
    return retval2;
}


int main()
{
    //return the results of each case
    int resolve = 0;
    std::cout << std::string(32,'-') << std::endl;
    resolve |= monolithicCase();
    std::cout << std::string(32,'-') << std::endl;
    resolve |= parserCase();
    std::cout << std::string(32,'-') << std::endl;
    resolve |= hardwiredCase();

    if (resolve)
    {
        std::cout << "Failure!\n";
    }
    else
    {
        std::cout << "Success!\n";
    }
    
    
    return resolve;
}

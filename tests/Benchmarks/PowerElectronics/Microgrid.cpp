

#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/Benchmark.hpp>

using index_type = size_t;
using real_type  = double;

/**
 * @brief Tests network of distributed generators.
 *
 * @param Nsize - The number of DG line load cobinations to generate for scale
 * @param error_tol - The tolerance for the model to meet to pass
 * @param debug_output - Enable debug output
 * @param use_DAE_keys - Choice between using DAE or ODE keys
 * @return int returns 0 if successful, >0 otherwise
 */
std::unique_ptr<GridKit::PowerElectronicsModel<double, size_t>> microgrid(index_type Nsize)
{
  using namespace GridKit;

  bool use_jac = true;

  real_type rel_tol = 1e-8;
  real_type abs_tol = 1e-8;

  // Create circuit model
  auto sys_model = std::make_unique<PowerElectronicsModel<real_type, index_type>>(rel_tol,
                                                                                  abs_tol,
                                                                                  true,
                                                                                  10000);

  // Modeled after the problem in the paper
  // Every Bus has the same virtual resistance. This is due to the numerical stability as mentioned in the paper.
  real_type RN = 1.0e4;

  // DG Params Vector
  // All DGs have the same set of parameters except for the first two.
  GridKit::DistributedGeneratorParameters<real_type, index_type> DG_parms1;
  DG_parms1.wb_  = 2.0 * M_PI * 50.0;
  DG_parms1.wc_  = 31.41;
  DG_parms1.mp_  = 9.4e-5;
  DG_parms1.Vn_  = 380.0;
  DG_parms1.nq_  = 1.3e-3;
  DG_parms1.F_   = 0.75;
  DG_parms1.Kiv_ = 420.0;
  DG_parms1.Kpv_ = 0.1;
  DG_parms1.Kic_ = 2.0e4;
  DG_parms1.Kpc_ = 15.0;
  DG_parms1.Cf_  = 5.0e-5;
  DG_parms1.rLf_ = 0.1;
  DG_parms1.Lf_  = 1.35e-3;
  DG_parms1.rLc_ = 0.03;
  DG_parms1.Lc_  = 0.35e-3;

  GridKit::DistributedGeneratorParameters<real_type, index_type> DG_parms2;
  DG_parms2.wb_  = 2.0 * M_PI * 50.0;
  DG_parms2.wc_  = 31.41;
  DG_parms2.mp_  = 12.5e-5;
  DG_parms2.Vn_  = 380.0;
  DG_parms2.nq_  = 1.5e-3;
  DG_parms2.F_   = 0.75;
  DG_parms2.Kiv_ = 390.0;
  DG_parms2.Kpv_ = 0.05;
  DG_parms2.Kic_ = 16.0e3;
  DG_parms2.Kpc_ = 10.5;
  DG_parms2.Cf_  = 50.0e-6;
  DG_parms2.rLf_ = 0.1;
  DG_parms2.Lf_  = 1.35e-3;
  DG_parms2.rLc_ = 0.03;
  DG_parms2.Lc_  = 0.35e-3;

  std::vector<GridKit::DistributedGeneratorParameters<real_type, index_type>> DGParams_list(2 * Nsize, DG_parms2);

  DGParams_list[0] = DG_parms1;
  DGParams_list[1] = DG_parms1;

  // line vector params
  // Every odd line has the same parameters and every even line has the same parameters
  real_type              rline1 = 0.23;
  real_type              Lline1 = 0.1 / (2.0 * M_PI * 50.0);
  real_type              rline2 = 0.35;
  real_type              Lline2 = 0.58 / (2.0 * M_PI * 50.0);
  std::vector<real_type> rline_list(2 * Nsize - 1, 0.0);
  std::vector<real_type> Lline_list(2 * Nsize - 1, 0.0);
  for (index_type i = 0; i < rline_list.size(); i++)
  {
    rline_list[i] = (i % 2) ? rline2 : rline1;
    Lline_list[i] = (i % 2) ? Lline2 : Lline1;
  }

  // load parms
  // Only the first load has the same paramaters.
  real_type rload1 = 3.0;
  real_type Lload1 = 2.0 / (2.0 * M_PI * 50.0);
  real_type rload2 = 2.0;
  real_type Lload2 = 1.0 / (2.0 * M_PI * 50.0);

  std::vector<real_type> rload_list(Nsize, rload2);
  std::vector<real_type> Lload_list(Nsize, Lload2);
  rload_list[0] = rload1;
  Lload_list[0] = Lload1;

  //							DGs	+		- refframe	   Lines +				Loads
  index_type vec_size_internals = 13 * (2 * Nsize) - 1 + (2 + 4 * (Nsize - 1)) + 2 * Nsize;
  //							\omegaref + BusDQ
  index_type vec_size_externals = 1 + 2 * (2 * Nsize);

  std::vector<index_type> vdqbus_index(2 * Nsize, 0);
  vdqbus_index[0] = vec_size_internals + 1;
  for (index_type i = 1; i < vdqbus_index.size(); i++)
  {
    vdqbus_index[i] = vdqbus_index[i - 1] + 2;
  }

  // Total size of the vector setup
  index_type vec_size_total = vec_size_internals + vec_size_externals;

  // Create the reference DG
  auto* dg_ref = new DistributedGenerator<real_type, index_type>(0,
                                                                 DGParams_list[0],
                                                                 true);
  // ref motor
  dg_ref->setExternalConnectionNodes(0, vec_size_internals);
  // outputs
  dg_ref->setExternalConnectionNodes(1, vdqbus_index[0]);
  dg_ref->setExternalConnectionNodes(2, vdqbus_index[0] + 1);
  //"grounding" of the difference
  dg_ref->setExternalConnectionNodes(3, static_cast<index_type>(-1));
  // internal connections
  for (index_type i = 0; i < 12; i++)
  {
    dg_ref->setExternalConnectionNodes(4 + i, i);
  }
  sys_model->addComponent(dg_ref);

  // Keep track of models and index location
  index_type indexv   = 12;
  index_type model_id = 1;
  // Add all other DGs
  for (index_type i = 1; i < 2 * Nsize; i++)
  {
    // current DG to add
    auto* dg = new DistributedGenerator<real_type, index_type>(model_id++,
                                                               DGParams_list[i],
                                                               false);
    // ref motor
    dg->setExternalConnectionNodes(0, vec_size_internals);
    // outputs
    dg->setExternalConnectionNodes(1, vdqbus_index[i]);
    dg->setExternalConnectionNodes(2, vdqbus_index[i] + 1);
    // internal connections
    for (index_type j = 0; j < 13; j++)
    {
      dg->setExternalConnectionNodes(3 + j, indexv + j);
    }
    indexv += 13;
    sys_model->addComponent(dg);
  }

  // Load all the Line compoenents
  for (index_type i = 0; i < 2 * Nsize - 1; i++)
  {
    // line
    auto* line_model = new MicrogridLine<real_type, index_type>(model_id++,
                                                                rline_list[i],
                                                                Lline_list[i]);
    // ref motor
    line_model->setExternalConnectionNodes(0, vec_size_internals);
    // input connections
    line_model->setExternalConnectionNodes(1, vdqbus_index[i]);
    line_model->setExternalConnectionNodes(2, vdqbus_index[i] + 1);
    // output connections
    line_model->setExternalConnectionNodes(3, vdqbus_index[i + 1]);
    line_model->setExternalConnectionNodes(4, vdqbus_index[i + 1] + 1);
    // internal connections
    for (index_type j = 0; j < 2; j++)
    {
      line_model->setExternalConnectionNodes(5 + j, indexv + j);
    }
    indexv += 2;
    sys_model->addComponent(line_model);
  }

  //  Load all the Load components
  for (index_type i = 0; i < Nsize; i++)
  {
    auto* load_model = new MicrogridLoad<real_type, index_type>(model_id++,
                                                                rload_list[i],
                                                                Lload_list[i]);
    // ref motor
    load_model->setExternalConnectionNodes(0, vec_size_internals);
    // input connections
    load_model->setExternalConnectionNodes(1, vdqbus_index[2 * i]);
    load_model->setExternalConnectionNodes(2, vdqbus_index[2 * i] + 1);
    // internal connections
    for (index_type j = 0; j < 2; j++)
    {
      load_model->setExternalConnectionNodes(3 + j, indexv + j);
    }
    indexv += 2;
    sys_model->addComponent(load_model);
  }

  // Add all the microgrid Virtual DQ Buses
  for (index_type i = 0; i < 2 * Nsize; i++)
  {
    auto* virDQbus_model = new MicrogridBusDQ<real_type, index_type>(model_id++, RN);

    virDQbus_model->setExternalConnectionNodes(0, vdqbus_index[i]);
    virDQbus_model->setExternalConnectionNodes(1, vdqbus_index[i] + 1);
    sys_model->addComponent(virDQbus_model);
  }

  // allocate all the intial conditions
  sys_model->allocate(vec_size_total);

  // Create Intial points for states. Every state is to specified to the zero intially
  for (index_type i = 0; i < vec_size_total; i++)
  {
    sys_model->y()[i]  = 0.0;
    sys_model->yp()[i] = 0.0;
  }

  // Create Intial derivatives specifics generated in MATLAB
  for (index_type i = 0; i < 2 * Nsize; i++)
  {
    sys_model->yp()[13 * i - 1 + 3] = DGParams_list[i].Vn_;
    sys_model->yp()[13 * i - 1 + 5] = DGParams_list[i].Kpv_ * DGParams_list[i].Vn_;
    sys_model->yp()[13 * i - 1 + 7] = (DGParams_list[i].Kpc_ * DGParams_list[i].Kpv_ * DGParams_list[i].Vn_) / DGParams_list[i].Lf_;
  }

  // since the intial P_com = 0, the set the intial vector to the reference frame
  sys_model->y()[vec_size_internals] = DG_parms1.wb_;

  sys_model->initialize();

  sys_model->updateTime(0.0, 1.0e-8);

  return std::move(sys_model);
}

void doIterations(index_type Nsize, bool use_csr)
{
  const size_t NUM_ITERATIONS = 20;

  real_type t_init  = 0.0;
  real_type t_final = 1.0;

  for (unsigned i = 0; i < NUM_ITERATIONS; i++)
  {
    std::cerr << i << ' ';
    GridKit::Utility::benchmark.newIteration();

    auto sys_model = microgrid(Nsize);
    std::cerr << std::format("Num components: {}, System Size: {}\n", sys_model->numComponents(), sys_model->size());

    // Create numerical integrator and configure it for the generator model
    auto* idas = new AnalysisManager::Sundials::Ida<real_type, index_type>(sys_model.get());
    idas->setUseCsr(use_csr);

    // setup simulation
    idas->configureSimulation();
    idas->getDefaultInitialCondition();
    idas->initializeSimulation(t_init);

    auto timer = GridKit::Utility::startTime<true>("Microgrid Simulation");
    idas->runSimulation(t_final);
    GridKit::Utility::endTime<true>(std::move(timer));
  }
  std::cerr << '\n';
}

/**
 * @brief Run Scale Microgrid test cases of N = (2,4,8) and check for correctness.
 *
 * @param argc unused
 * @param argv unsued
 * @return int
 */
int main(int /* argc */, char const** /* argv */)
{
  GridKit::Utility::benchmark.newRun("Microgrid N=2 COO");
  doIterations(2, false);

  GridKit::Utility::benchmark.newRun("Microgrid N=2 CSR");
  doIterations(2, true);

  GridKit::Utility::benchmark.newRun("Microgrid N=64 COO");
  doIterations(64, false);

  GridKit::Utility::benchmark.newRun("Microgrid N=64 CSR");
  doIterations(64, true);

  std::cerr << GridKit::Utility::benchmark.report() << '\n';
  return 0;
}
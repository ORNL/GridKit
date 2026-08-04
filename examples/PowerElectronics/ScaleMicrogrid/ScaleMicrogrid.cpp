#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/SignalNode.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

using index_type = size_t;
using real_type  = double;

// Include solution keys for the three test cases N = (2, 4, 8) plus tolerances here:
#include "SolutionKeys.hpp"

static int test(index_type Nsize, real_type test_tolerance, bool error_tol = false);

/**
 * @brief Run Scale Microgrid test cases of N = (2,4,8) and check for correctness.
 *
 * @param argc unused
 * @param argv unsued
 * @return int
 */
int main(int /* argc */, char const** /* argv */)
{
  int       retval    = 0;
  bool      debug_out = false;
  real_type tol       = SCALE_MICROGRID_ERROR_TOL;

  retval += test(2, tol, debug_out);
  retval += test(4, tol, debug_out);
  retval += test(8, tol, debug_out);
  if (retval > 0)
  {
    std::cout << "Some tests fail!!\n";
  }
  else
  {
    std::cout << "All tests pass!!\n";
  }
  return retval;
}

/**
 * @brief Tests network of distributed generators.
 *
 * @param Nsize - The number of DG line load cobinations to generate for scale
 * @param error_tol - The tolerance for the model to meet to pass
 * @param debug_output - Enable debug output
 * @param use_DAE_keys - Choice between using DAE or ODE keys
 * @return int returns 0 if successful, >0 otherwise
 */
int test(index_type Nsize, real_type error_tol, bool debug_output)
{
  using namespace GridKit;

  bool use_jac = true;

  real_type t_init  = 0.0;
  real_type t_final = 1.0;

  real_type rel_tol = SCALE_MICROGRID_REL_TOL;
  real_type abs_tol = SCALE_MICROGRID_ABS_TOL;

  // Create circuit model
  auto* sys_model = new PowerElectronicsModel<real_type, index_type>(use_jac);

  const std::vector<real_type>* true_vec = &answer_key_N8;

  switch (Nsize)
  {
  case 2:
    true_vec = &answer_key_N2;
    break;
  case 4:
    true_vec = &answer_key_N4;
    break;
  case 8:
    true_vec = &answer_key_N8;
    break;
  default:
    std::cout << "No reference solution for Nsize = " << Nsize << ".\n";
    std::cout << "Using default Nsize = 8.\n";
  }

  // Modeled after the problem in the paper
  // Every Bus has the same virtual resistance. This is due to the numerical stability as mentioned in the paper.
  real_type RN = 1.0e4;

  // DG Params Vector
  // All DGs have the same set of parameters except for the first two.
  static constexpr auto pi = std::numbers::pi_v<real_type>;
  GridKit::DistributedGeneratorParameters<real_type, index_type> DG_parms1;
  DG_parms1.wb_  = 2.0 * pi * 50.0;
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
  DG_parms2.wb_  = 2.0 * pi * 50.0;
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
  real_type              Lline1 = 0.1 / (2.0 * pi * 50.0);
  real_type              rline2 = 0.35;
  real_type              Lline2 = 0.58 / (2.0 * pi * 50.0);
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
  real_type Lload1 = 2.0 / (2.0 * pi * 50.0);
  real_type rload2 = 2.0;
  real_type Lload2 = 1.0 / (2.0 * pi * 50.0);

  std::vector<real_type> rload_list(Nsize, rload2);
  std::vector<real_type> Lload_list(Nsize, Lload2);
  rload_list[0] = rload1;
  Lload_list[0] = Lload1;

  using SignalNode = GridKit::PowerElectronics::SignalNode<double, size_t>;
  SignalNode dg_signal;
  sys_model->addNode(&dg_signal);

  using Bus                                     = GridKit::PowerElectronics::MicrogridBus<double, size_t>;
  std::unique_ptr<std::unique_ptr<Bus>[]> buses = std::make_unique<std::unique_ptr<Bus>[]>(2 * Nsize);
  for (size_t i = 0; i < 2 * Nsize; i++)
  {
    buses[i] = std::make_unique<Bus>();
    sys_model->addNode(buses[i].get());
  }

  // Create the reference DG
  auto* dg_ref = new DistributedGenerator<real_type, index_type>(0,
                                                                 DGParams_list[0],
                                                                 true,
                                                                 &dg_signal,
                                                                 buses[0].get());
  sys_model->addComponent(dg_ref);

  // Keep track of models and index location
  index_type model_id = 1;
  // Add all other DGs
  for (index_type i = 1; i < 2 * Nsize; i++)
  {
    // current DG to add
    auto* dg = new DistributedGenerator<real_type, index_type>(model_id++,
                                                               DGParams_list[i],
                                                               false,
                                                               &dg_signal,
                                                               buses[i].get());
    sys_model->addComponent(dg);
  }

  // Load all the Line compoenents
  for (index_type i = 0; i < 2 * Nsize - 1; i++)
  {
    // line
    auto* line_model = new MicrogridLine<real_type, index_type>(model_id++,
                                                                rline_list[i],
                                                                Lline_list[i],
                                                                &dg_signal,
                                                                buses[i].get(),
                                                                buses[i + 1].get());
    sys_model->addComponent(line_model);
  }

  //  Load all the Load components
  for (index_type i = 0; i < Nsize; i++)
  {
    auto* load_model = new MicrogridLoad<real_type, index_type>(model_id++,
                                                                rload_list[i],
                                                                Lload_list[i],
                                                                &dg_signal,
                                                                buses[2 * i].get());
    sys_model->addComponent(load_model);
  }

  // Add all the microgrid Virtual DQ Buses
  for (index_type i = 0; i < 2 * Nsize; i++)
  {
    auto* virDQbus_model = new MicrogridBusDQ<real_type, index_type>(model_id++, RN, buses[i].get());
    sys_model->addComponent(virDQbus_model);
  }

  // allocate all the intial conditions
  sys_model->allocate();

  if (debug_output)
  {
    std::cout << sys_model->y().getSize() << std::endl;
  }

  auto* y  = sys_model->y().getData();
  auto* yp = sys_model->yp().getData();

  // Create initial points for states. Every state is to specified to the zero intially
  for (index_type i = 0; i < sys_model->size(); i++)
  {
    y[i]  = 0.0;
    yp[i] = 0.0;
  }

  // Create initial derivatives specifics generated in MATLAB
  for (index_type i = 0; i < 2 * Nsize; i++)
  {
    yp[13 * i - 1 + 2] = DGParams_list[i].Vn_;
    yp[13 * i - 1 + 4] = DGParams_list[i].Kpv_ * DGParams_list[i].Vn_;
    yp[13 * i - 1 + 6] = (DGParams_list[i].Kpc_ * DGParams_list[i].Kpv_ * DGParams_list[i].Vn_) / DGParams_list[i].Lf_;
  }

  // since the intial P_com = 0, the set the intial vector to the reference frame
  y[dg_signal.getNodeConnection(0).idx_] = DG_parms1.wb_;

  sys_model->y().setDataUpdated();
  sys_model->yp().setDataUpdated();

  sys_model->initialize();
  sys_model->evaluateResidual();
  auto&       fres      = sys_model->getResidual();
  const auto* fres_data = fres.getData();
  if (debug_output)
  {
    std::cout << "Verify initial resisdual is zero: {\n";
    for (index_type i = 0; i < fres.getSize(); i++)
    {
      std::cout << i << " : " << fres_data[i] << "\n";
    }
    std::cout << "}\n";
  }

  sys_model->updateTime(0.0, 1.0e-8);
  sys_model->evaluateJacobian();

  if (debug_output)
  {
    std::cout << "Initial Jacobian with alpha:\n";
    sys_model->getCsrJacobian()->print();
  }

  // Create numerical integrator and configure it for the generator model
  auto* idas = new AnalysisManager::Sundials::Ida<real_type, index_type>(sys_model);

  // setup simulation
  idas->setTolerance(rel_tol, abs_tol);
  idas->setMaxSteps(SCALE_MICROGRID_MAX_STEPS);
  idas->configureSimulation();
  idas->getDefaultInitialCondition();
  idas->initializeSimulation(t_init);

  idas->runSimulation(t_final);

  auto&       yfinal      = sys_model->y();
  const auto* yfinal_data = yfinal.getData();

  if (debug_output)
  {
    std::cout << "Final Vector y\n";
    for (index_type i = 0; i < yfinal.getSize(); i++)
    {
      std::cout << i << " : " << yfinal_data[i] << "\n";
    }
  }

  bool test_pass = true;

  real_type sum_top    = 0.0;
  real_type sum_bottom = 0.0;

  // check relative error
  std::cout << "Test the relative error for N = " << Nsize << "\n";
  for (index_type i = 0; i < true_vec->size(); i++)
  {
    // Print the Elementwise Relative Error
    if (debug_output)
      std::cout << i << " : " << abs(true_vec->at(i) - yfinal_data[i]) / abs(true_vec->at(i)) << "\n";

    sum_top    += (true_vec->at(i) - yfinal_data[i]) * (true_vec->at(i) - yfinal_data[i]);
    sum_bottom += (true_vec->at(i) * true_vec->at(i));
  }

  real_type norm2error = (std::sqrt(sum_top) / std::sqrt(sum_bottom));
  std::cout << "2-Norm relative error: " << norm2error << std::endl;
  test_pass = norm2error < error_tol;

  delete idas;
  delete sys_model;

  if (test_pass)
  {
    std::cout << "Test with Nsize = " << Nsize << " passes!\n";
    return 0;
  }
  else
  {
    std::cout << "Test with Nsize = " << Nsize << " fails!\n";
    return 1;
  }
}

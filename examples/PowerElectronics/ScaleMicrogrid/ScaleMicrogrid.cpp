#include <cmath>
#include <iostream>
#include <vector>

#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

#include <examples/PowerElectronics/PowerElectronicsExamplesHelper/MicrogridNetwork.hpp>

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

  // Build and assemble the scaled microgrid network.
  ScaleMicrogridNetwork<double, size_t> network(Nsize);
  assembleSystem(network, *sys_model);

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
    const auto& params = network.DGParam_list[i];

    yp[13 * i - 1 + 2] = params.Vn_;
    yp[13 * i - 1 + 4] = params.Kpv_ * params.Vn_;
    yp[13 * i - 1 + 6] = (params.Kpc_ * params.Kpv_ * params.Vn_) / params.Lf_;
  }

  // since the intial P_com = 0, the set the intial vector to the reference frame
  y[network.dg_signal.getNodeConnection(0).idx_] = network.DGParam_list[0].wb_;

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

#include <cmath>
#include <iostream>

#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "Common/MicrogridNetwork.hpp"

using index_type = size_t;
using real_type  = double;

// Include solution keys for the three test cases N = (2, 4, 8) plus tolerances here:
#include "SolutionKeys.hpp"

int printMicrogridSystems(index_type N_size);

/**
 * @brief Run Scale Microgrid for a specific N given by the user.
 *
 * @param[in] argc should be 1 if no N given, 2 if N given
 * @param[in] argv can contain the N value (argv[1])
 * @return int
 */
int main(int argc, char const* argv[])
{
  // Default value
  index_type N_size = 4;

  // Parse command line arguments if provided
  if (argc > 1)
  {
    try
    {
      N_size = static_cast<index_type>(std::stoi(argv[1]));
    }
    catch (const std::exception& e)
    {
      std::cerr << "Error parsing grid size argument: " << e.what() << std::endl;
      std::cerr << "Using default grid size = " << N_size << std::endl;
    }
  }

  return printMicrogridSystems(N_size);
}

/**
 * @brief Tests network of distributed generators.
 *
 * @param[in] N_size - The number of DG line load cobinations to generate for scale
 * @return int returns 0 if successful, >0 otherwise
 */
int printMicrogridSystems(index_type N_size)
{
  using namespace GridKit;

  bool use_jac = true;

  real_type t_init  = 0.0;
  real_type t_final = 1.0;

  real_type rel_tol = SCALE_MICROGRID_REL_TOL;

  // Create circuit model
  PowerElectronicsModel<real_type, index_type> sys_model(use_jac);

  // Ensure minimum size requirement
  if (N_size < 1)
  {
    std::cout << "N_size must be at least 1.\n";
    return 1;
  }

  // Build and assemble the scaled microgrid network.
  ScaleMicrogridNetwork network(N_size);

  buildScaleMicrogridNetwork(network);
  assembleSystem(network, sys_model);

  // allocate all the initial conditions
  sys_model.allocate();

  auto* y  = sys_model.y().getData();
  auto* yp = sys_model.yp().getData();

  // Create Initial points for states. Every state is set to zero initially
  for (index_type i = 0; i < sys_model.size(); i++)
  {
    y[i]  = 0.0;
    yp[i] = 0.0;
  }

  //-------------------------------------------------------------------
  // Create Initial derivatives specifics generated in MATLAB
  //-------------------------------------------------------------------
  for (index_type i = 0; i < 2 * N_size; i++)
  {
    const auto& params = network.DGParam_list[i];

    yp[13 * i - 1 + 3] = params.Vn_;
    yp[13 * i - 1 + 5] = params.Kpv_ * params.Vn_;
    yp[13 * i - 1 + 7] = (params.Kpc_ * params.Kpv_ * params.Vn_) / params.Lf_;
  }

  //---------------------------------------------------------------------------
  // since the initial P_com = 0, set the initial vector to the reference frame
  //---------------------------------------------------------------------------
  y[network.dg_signal.getNodeConnection(0).idx_] = network.DGParam_list[0].wb_;

  sys_model.y().setDataUpdated();
  sys_model.yp().setDataUpdated();

  sys_model.initialize();
  sys_model.evaluateResidual();

  // Output file names based on grid size
  std::string size_suffix = std::to_string(N_size);

  sys_model.updateTime(0.0, 1.0e-8);
  sys_model.evaluateJacobian();

  // Create numerical integrator and configure it for the generator model
  AnalysisManager::Sundials::Ida<real_type, index_type> idas(&sys_model);

  // setup simulation
  idas.setTolerance(rel_tol);
  idas.setMaxSteps(SCALE_MICROGRID_MAX_STEPS);
  idas.configureSimulation();
  idas.getDefaultInitialCondition();
  idas.initializeSimulation(t_init);
  idas.runSimulation(t_final);

  return 0;
}
#include <cmath>
#include <cstddef>
#include <iostream>

#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include <examples/PowerElectronics/PowerElectronicsExamplesHelper/MicrogridNetwork.hpp>

int main(int /* argc */, char const** /* argv */)
{
  /// @todo Needs to be modified. Some components are small relative to others thus
  /// there error is high (or could be matlab vector issue)
  double rel_tol         = 1.0e-8;
  size_t max_step_number = 3000;
  bool   use_jac         = true;
  bool   debug_output    = true;

  // Create model
  auto* sysmodel = new GridKit::PowerElectronicsModel<double, size_t>(use_jac);

  // Build the four-generator microgrid network.
  size_t                                N_size = 2;
  ScaleMicrogridNetwork<double, size_t> network(N_size);
  assembleSystem(network, *sysmodel);

  // Generator parameters used to construct the initial conditions.
  const auto& parms1 = network.DGParam_list[0];
  const auto& parms2 = network.DGParam_list[2];

  sysmodel->allocate();

  std::cout << sysmodel->y().getSize() << std::endl;

  auto* y  = sysmodel->y().getData();
  auto* yp = sysmodel->yp().getData();

  // Create initial points for states
  for (size_t i = 0; i < sysmodel->size(); i++)
  {
    y[i]  = 0.0;
    yp[i] = 0.0;
  }

  // Create initial derivatives specifics generated in MATLAB
  // DGs 1
  yp[2]      = parms1.Vn_;
  yp[4]      = parms1.Kpv_ * parms1.Vn_;
  yp[6]      = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
  yp[12 + 2] = parms1.Vn_;
  yp[12 + 4] = parms1.Kpv_ * parms1.Vn_;
  yp[12 + 6] = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
  for (size_t i = 2; i < 4; i++)
  {
    yp[13 * i - 1 + 2] = parms2.Vn_;
    yp[13 * i - 1 + 4] = parms2.Kpv_ * parms2.Vn_;
    yp[13 * i - 1 + 6] = (parms2.Kpc_ * parms2.Kpv_ * parms2.Vn_) / parms2.Lf_;
  }

  // since the intial P_com = 0
  y[network.dg_signal.getNodeConnection(0).idx_] = parms1.wb_;

  sysmodel->y().setDataUpdated();
  sysmodel->yp().setDataUpdated();

  sysmodel->initialize();
  sysmodel->evaluateResidual();

  // Optional debuging output
  if (debug_output)
  {
    auto&       fres      = sysmodel->getResidual();
    const auto* fres_data = fres.getData();
    std::cout << "Verify initial resisdual is zero: {\n";
    for (size_t i = 0; i < fres.getSize(); i++)
    {
      std::cout << i << " :" << fres_data[i] << "\n";
    }
    std::cout << "}\n";
  }

  sysmodel->updateTime(0.0, 1.0e-8);
  sysmodel->evaluateJacobian();

  // Optional debuging output
  if (debug_output)
  {
    std::cout << "Initial Jacobian with alpha:\n";
    sysmodel->getCsrJacobian()->print();
  }

  sysmodel->tagDifferentiable();

  bool all_internal_diff = true;
  bool all_external_alg  = true;

  const size_t num_node_vars = network.buses[0].size() + network.buses[1].size() + network.buses[2].size() + network.buses[3].size() + network.dg_signal.size();

  for (size_t i = 0; i < sysmodel->size() - num_node_vars; i++)
  {
    all_internal_diff = all_internal_diff && sysmodel->tag()[i];
    if (!sysmodel->tag()[i])
    {
      std::cout << "Unexepected algebraic-tagged internal variable found in index " << i << '\n';
    }
  }

  for (size_t i = sysmodel->size() - num_node_vars; i < sysmodel->size(); i++)
  {
    all_external_alg = all_external_alg && !sysmodel->tag()[i];
    if (sysmodel->tag()[i])
    {
      std::cout << "Unexepected differential-tagged external variable found in index " << i << '\n';
    }
  }

  std::cout << "Verify all internal variables are differential: " << all_internal_diff << ", and all external variables are algebraic: " << all_external_alg << '\n';

  // Create numerical integrator and configure it for the generator model
  auto* idas = new AnalysisManager::Sundials::Ida<double, size_t>(sysmodel);

  double t_init  = 0.0;
  double t_final = 1.0;

  // setup simulation
  idas->setTolerance(rel_tol);
  idas->setMaxSteps(max_step_number);
  idas->configureSimulation();
  idas->getDefaultInitialCondition();
  idas->initializeSimulation(t_init);

  idas->runSimulation(t_final);

  auto&       yfinal      = sysmodel->y();
  const auto* yfinal_data = yfinal.getData();

  // Optional debugging output
  if (debug_output)
  {
    std::cout << "Final vector y\n";
    for (size_t i = 0; i < yfinal.getSize(); i++)
    {
      std::cout << yfinal_data[i] << "\n";
    }
  }

  // Generated from MATLAB code ODE form with tolerances of 1e-12
  std::vector<double> true_vec{
      2.297543153595780e+04,
      1.275311524125022e+04,
      3.763060183116022e-02,
      -2.098153459325261e-02,
      1.848285659119097e-02,
      -1.563291404944864e-04,
      6.321941907011718e+01,
      -2.942264300846256e+01,
      3.634209302905854e+02,
      -2.668928293656362e-06,
      6.321941919221522e+01,
      -3.509200178595996e+01,
      2.297580486511343e+04,
      8.742028429066131e+03,
      3.710079564796484e-02,
      -1.421122598056797e-02,
      1.874079517807597e-02,
      -9.891304812687215e-05,
      6.232933298360234e+01,
      -1.796494061423331e+01,
      3.686353885026506e+02,
      3.465673854181523e-05,
      6.232933406188410e+01,
      -2.371564475187742e+01,
      -7.555954467454730e-03,
      1.727775042678524e+04,
      1.649365247247288e+04,
      3.116555157570849e-02,
      -2.985990066758010e-02,
      2.250012115906506e-02,
      -2.643873146501096e-04,
      4.861823510250247e+01,
      -4.088592755441309e+01,
      3.552597163751238e+02,
      -1.496407194199739e-04,
      4.861823504694532e+01,
      -4.642797132602495e+01,
      -8.273939686941580e-02,
      1.727723725566433e+04,
      9.182386962936238e+03,
      3.024959333190777e-02,
      -1.617250828202081e-02,
      2.318056864131751e-02,
      -1.295918667730514e-04,
      4.718938244522050e+01,
      -1.935782085675469e+01,
      3.662262287803608e+02,
      1.076423957830039e-04,
      4.718938116520511e+01,
      -2.507094256286497e+01,
      -8.445727984408551e-02,
      -1.881248349415025e+01,
      2.114714832305742e+01,
      4.329946674909793e+01,
      -3.037887936225145e+00,
      -4.487023117352992e+01,
      2.895883729832657e+01,
      8.199613345691378e+01,
      -5.623856502948122e+01,
      1.327498499660322e+02,
      -8.228065162347022e+01,
      3.119995747945993e+02,
      3.576922945168803e+02,
      -5.850795361581618e+00,
      3.641193316268954e+02,
      -8.846325267612976e+00,
      3.472146752739036e+02,
      -3.272400970143252e+01,
      3.604108939430972e+02,
      -3.492842627398574e+01};

  std::cout << "Testing Migrogrid ...\n";
  double error_allowed = 1e-4;
  double max_error     = 0.0;
  for (size_t i = 0; i < true_vec.size(); i++)
  {
    double error = std::abs(true_vec[i] - yfinal_data[i]) / std::abs(1.0 + true_vec[i]);
    if (error > max_error)
      max_error = error;
    if (error > error_allowed)
    {
      std::cout << "Model error for equation " << i << " is: " << error << "\n";
      std::cout << "Maximum allowed error is: " << error_allowed << "\n";
      std::cout << "Test FAILED!\n";
      return 1;
    }
  }
  std::cout << "Max error     = " << max_error << "\n";
  std::cout << "Allowed error = " << error_allowed << "\n";
  std::cout << "Test successful!\n";

  delete idas;
  delete sysmodel;

  return 0;
}

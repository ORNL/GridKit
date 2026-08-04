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

  // Modeled after the problem in the paper
  double RN = 1.0e4;

  // DG Params
  static constexpr auto pi = std::numbers::pi_v<double>;
  GridKit::DistributedGeneratorParameters<double, size_t> parms1;
  parms1.wb_  = 2.0 * pi * 50.0;
  parms1.wc_  = 31.41;
  parms1.mp_  = 9.4e-5;
  parms1.Vn_  = 380.0;
  parms1.nq_  = 1.3e-3;
  parms1.F_   = 0.75;
  parms1.Kiv_ = 420.0;
  parms1.Kpv_ = 0.1;
  parms1.Kic_ = 2.0e4;
  parms1.Kpc_ = 15.0;
  parms1.Cf_  = 5.0e-5;
  parms1.rLf_ = 0.1;
  parms1.Lf_  = 1.35e-3;
  parms1.rLc_ = 0.03;
  parms1.Lc_  = 0.35e-3;

  GridKit::DistributedGeneratorParameters<double, size_t> parms2;
  // Parameters from MATLAB Microgrid code for first DG
  parms2.wb_  = 2.0 * pi * 50.0;
  parms2.wc_  = 31.41;
  parms2.mp_  = 12.5e-5;
  parms2.Vn_  = 380.0;
  parms2.nq_  = 1.5e-3;
  parms2.F_   = 0.75;
  parms2.Kiv_ = 390.0;
  parms2.Kpv_ = 0.05;
  parms2.Kic_ = 16.0e3;
  parms2.Kpc_ = 10.5;
  parms2.Cf_  = 50.0e-6;
  parms2.rLf_ = 0.1;
  parms2.Lf_  = 1.35e-3;
  parms2.rLc_ = 0.03;
  parms2.Lc_  = 0.35e-3;

  // Line params
  double rline1 = 0.23;
  double Lline1 = 0.1 / (2.0 * pi * 50.0);

  double rline2 = 0.35;
  double Lline2 = 0.58 / (2.0 * pi * 50.0);

  double rline3 = 0.23;
  double Lline3 = 0.1 / (2.0 * pi * 50.0);

  // load parms
  double rload1 = 3.0;
  double Lload1 = 2.0 / (2.0 * pi * 50.0);

  double rload2 = 2.0;
  double Lload2 = 1.0 / (2.0 * pi * 50.0);

  using SignalNode = GridKit::PowerElectronics::SignalNode<double, size_t>;
  SignalNode dg_signal;

  sysmodel->addNode(&dg_signal);

  using Bus = GridKit::PowerElectronics::MicrogridBus<double, size_t>;
  Bus bus1;
  Bus bus2;
  Bus bus3;
  Bus bus4;

  sysmodel->addNode(&bus1);
  sysmodel->addNode(&bus2);
  sysmodel->addNode(&bus3);
  sysmodel->addNode(&bus4);

  // dg 1
  GridKit::DistributedGenerator<double, size_t>* dg1 = new GridKit::DistributedGenerator<double, size_t>(
      0, parms1, true, &dg_signal, &bus1);
  sysmodel->addComponent(dg1);

  // dg 2
  GridKit::DistributedGenerator<double, size_t>* dg2 = new GridKit::DistributedGenerator<double, size_t>(
      1, parms1, false, &dg_signal, &bus2);
  sysmodel->addComponent(dg2);

  // dg 3
  GridKit::DistributedGenerator<double, size_t>* dg3 = new GridKit::DistributedGenerator<double, size_t>(
      2, parms2, false, &dg_signal, &bus3);
  sysmodel->addComponent(dg3);

  // dg 4
  GridKit::DistributedGenerator<double, size_t>* dg4 = new GridKit::DistributedGenerator<double, size_t>(
      3, parms2, false, &dg_signal, &bus4);
  sysmodel->addComponent(dg4);

  // Lines

  // line 1
  GridKit::MicrogridLine<double, size_t>* l1 = new GridKit::MicrogridLine<double, size_t>(
      4, rline1, Lline1, &dg_signal, &bus1, &bus2);
  sysmodel->addComponent(l1);

  // line 2
  GridKit::MicrogridLine<double, size_t>* l2 = new GridKit::MicrogridLine<double, size_t>(
      5, rline2, Lline2, &dg_signal, &bus2, &bus3);
  sysmodel->addComponent(l2);

  // line 3
  GridKit::MicrogridLine<double, size_t>* l3 = new GridKit::MicrogridLine<double, size_t>(
      6, rline3, Lline3, &dg_signal, &bus3, &bus4);
  sysmodel->addComponent(l3);

  //  loads

  // load 1
  GridKit::MicrogridLoad<double, size_t>* load1 = new GridKit::MicrogridLoad<double, size_t>(7, rload1, Lload1, &dg_signal, &bus1);
  sysmodel->addComponent(load1);

  // load 2
  GridKit::MicrogridLoad<double, size_t>* load2 = new GridKit::MicrogridLoad<double, size_t>(8, rload2, Lload2, &dg_signal, &bus3);
  sysmodel->addComponent(load2);

  // Virtual PQ Buses
  GridKit::MicrogridBusDQ<double, size_t>* bus_para_1 = new GridKit::MicrogridBusDQ<double, size_t>(9, RN, &bus1);
  sysmodel->addComponent(bus_para_1);

  GridKit::MicrogridBusDQ<double, size_t>* bus_para_2 = new GridKit::MicrogridBusDQ<double, size_t>(10, RN, &bus2);
  sysmodel->addComponent(bus_para_2);

  GridKit::MicrogridBusDQ<double, size_t>* bus_para_3 = new GridKit::MicrogridBusDQ<double, size_t>(11, RN, &bus3);
  sysmodel->addComponent(bus_para_3);

  GridKit::MicrogridBusDQ<double, size_t>* bus_para_4 = new GridKit::MicrogridBusDQ<double, size_t>(12, RN, &bus4);
  sysmodel->addComponent(bus_para_4);

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
  y[dg_signal.getNodeConnection(0).idx_] = parms1.wb_;

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

/**
 * Order test of MultiStageCosimulation on the 4-variable partitioned DAE.
 *
 *   x' = y - x^3 + z
 *   y' = -x + y(1 - y^2) - w
 *   0  = sin(x - y) - z
 *   0  = w + y - z - pi/2
 *
 * 1. Solve the DAE monolithically with tight IDA tolerances (reference), and
 *    sanity-check it against the known solution at T_END.
 * 2. Solve it with MultiStageCosimulation using N = 10, 20, 40, ... macro steps.
 * 3. Print the max error at T_END and the observed order p = log2(e_N / e_2N).
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Dynamic/MultiStageCosimulation.hpp>

#include "tests/UnitTests/PowerElectronics/CosimPartitionedDAETests.hpp"

using Bus       = GridKit::PowerElectronics::MicrogridBus<double, int>;
using System    = GridKit::PowerElectronicsModel<double, int>;
using Subsystem = GridKit::SubsystemModel<double, int>;
using Interface = GridKit::PowerElectronics::BusPartitionInterface<double, int>;
using Cosim     = AnalysisManager::MultiStageCosimulation<double, int>;
using Ida       = AnalysisManager::Sundials::Ida<double, int>;

const double T_END = 5.0;
const double PI    = GridKit::PartitionedDae::PI;

// Global ordering: component internals first, then bus node variables.
enum GlobalIndex
{
  IX = 0, // x (component, differential)
  IZ = 1, // z (component, algebraic)
  IY = 2, // y (bus, differential)
  IW = 3  // w (bus, algebraic)
};

/**
 * Partitioned DAE assembled as one system. Global ordering: x z y w.
 */
struct PartitionedDae
{
  Bus                                            bus;
  System*                                        system;
  GridKit::PartitionedDaeComponent<double, int>* comp;
  GridKit::PartitionedDaeBus<double, int>*       bus_comp;

  PartitionedDae()
  {
    system   = new System();
    comp     = new GridKit::PartitionedDaeComponent<double, int>(&bus, 1);
    bus_comp = new GridKit::PartitionedDaeBus<double, int>(&bus, 2);

    system->addComponent(comp);
    system->addComponent(bus_comp);
    system->addNode(&bus);
    system->allocate();

    // bus.tag()[0] = true;

    // Consistent initial state: (x, y, z, w) = (pi, pi/2, 1, 1)
    double* y  = system->y().getData();
    double* yp = system->yp().getData();
    for (int i = 0; i < system->size(); ++i)
    {
      y[i]  = 0.0;
      yp[i] = 0.0;
    }
    y[IX] = PI;
    y[IZ] = 1.0;
    y[IY] = PI / 2.0;
    y[IW] = 1.0;

    // Differential rows: with y' = 0 the residual F = g(y) - y' equals g(y).
    system->updateTime(0.0, 0.0);
    system->evaluateResidual();
    const double* F = system->getResidual().getData();
    yp[IX]          = F[IX];
    yp[IY]          = F[IY];

    // Algebraic rows: differentiate the constraints so y' is fully consistent
    //   z = sin(x - y)       => z' = cos(x - y) (x' - y')
    //   w = z - y + pi/2     => w' = z' - y'
    yp[IZ] = std::cos(y[IX] - y[IY]) * (yp[IX] - yp[IY]);
    yp[IW] = yp[IZ] - yp[IY];

    system->y().setDataUpdated();
    system->yp().setDataUpdated();
  }

  ~PartitionedDae()
  {
    delete system;
  }
};

std::vector<double> solveMonolithic()
{
  PartitionedDae dae;

  Ida ida(dae.system);
  ida.setTolerance(1e-14, 1e-14);
  ida.setMaxSteps(1000000);
  ida.configureSimulation();
  ida.getDefaultInitialCondition();
  ida.initializeSimulation(0.0, false);
  ida.runSimulation(T_END);

  const double* y = dae.system->y().getData();
  return std::vector<double>(y, y + dae.system->size());
}

std::vector<double> solveCosim(int n_steps)
{
  PartitionedDae dae;

  // Partition 1: Bus (+ interface carrying the component's bus terms)
  // Partition 2: Component (x, z)
  auto* bus_interface = new Interface(&dae.bus, dae.comp, 3);
  bus_interface->allocate();

  auto* partition1 = new Subsystem();
  partition1->addComponent(dae.bus_comp);
  partition1->addInterface(bus_interface);
  partition1->addNode(&dae.bus);

  auto* partition2 = new Subsystem();
  partition2->addComponent(dae.comp);

  // Macro time points 0, dt, 2dt, ..., T_END
  std::vector<double> time_points(n_steps + 1);
  for (int k = 0; k <= n_steps; ++k)
  {
    time_points[k] = T_END * k / n_steps;
  }

  size_t num_iterations = 1;
  size_t num_stages     = 4;
  double  tol           = 1e-12;

  auto* cosim = new Cosim(num_iterations, num_stages, tol);
  cosim->addSubsystems(partition1);
  cosim->addSubsystems(partition2);
  cosim->allocate();
  cosim->initializeSimulation(time_points, &dae.system->y(), &dae.system->yp());
  cosim->runSimulation();

  const double*       d = cosim->y().getData(GridKit::memory::HOST);
  std::vector<double> y(d, d + cosim->y().getSize());

  // cosim is not deleted: its inline destructor needs ReSolve at link time
  delete partition1;
  delete partition2;
  return y;
}

int main()
{
  const std::vector<double> ref = solveMonolithic();

  std::cout << "Monolithic system solved!" << std::endl;
  // Sanity check: monolithic solution vs. known solution at T_END
  std::vector<double> known(4);
  known[IX] = -1.017881864730135e+00;
  known[IY] = -1.641238691143580e+00;
  known[IZ] = 5.837639304756466e-01;
  known[IW] = 3.795798948414124e+00;

  double ref_error = 0.0;
  for (size_t i = 0; i < known.size(); ++i)
  {
    ref_error = std::max(ref_error, std::abs(ref[i] - known[i]) / (std::abs(known[i]) + 1));
  }

  const std::vector<int> step_counts = {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5000};
  std::vector<double>    errors;

  for (int n : step_counts)
  {
    const std::vector<double> y = solveCosim(n);

    double error = 0.0;
    for (size_t i = 0; i < ref.size(); ++i)
    {
      error = std::max(error, std::abs(y[i] - ref[i]) / (std::abs(ref[i]) + 1));
    }
    errors.push_back(error);
  }

  // Printed at the end so it is not buried in the cosimulation's own output
  std::cout << "\nMonolithic reference vs. known solution at T_END: "
            << std::scientific << std::setprecision(3) << ref_error << "\n";

  std::cout << "\n    N          dt     max error     p\n";
  for (size_t k = 0; k < step_counts.size(); ++k)
  {
    std::cout << std::setw(5) << step_counts[k]
              << std::scientific << std::setprecision(3)
              << std::setw(12) << T_END / step_counts[k]
              << std::setw(14) << errors[k];

    if (k > 0)
    {
      std::cout << std::fixed << std::setprecision(2)
                << std::setw(6) << std::log2(errors[k - 1] / errors[k]);
    }
    std::cout << "\n";
  }

  std::cout << "\n        dt     max error\n";
  for (size_t k = 0; k < step_counts.size(); ++k)
  {
    std::cout << std::scientific << std::setprecision(3)
              << std::setw(12) << T_END / step_counts[k]
              << std::setw(14) << errors[k]
              << std::endl;
  }

  return 0;
}
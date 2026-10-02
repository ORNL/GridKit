/**
 * Order test of MultiStageCosimulation on the HIRES problem.
 *
 * 1. Solve HIRES monolithically with tight IDA tolerances (reference).
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

#include "CosimHiresTests.hpp"

using Bus       = GridKit::PowerElectronics::MicrogridBus<double, int>;
using System    = GridKit::PowerElectronicsModel<double, int>;
using Subsystem = GridKit::SubsystemModel<double, int>;
using Interface = GridKit::PowerElectronics::BusPartitionInterface<double, int>;
using Cosim     = AnalysisManager::MultiStageCosimulation<double, int>;
using Ida       = AnalysisManager::Sundials::Ida<double, int>;

const double T_END = 5.0;

/**
 * HIRES assembled as one system. Global ordering: y1 y2 y3 y6 y7 y8 y4 y5.
 */
struct Hires
{
  Bus                                    bus;
  System*                                system;
  GridKit::HiresComponent1<double, int>* comp1;
  GridKit::HiresBus<double, int>*        bus_comp;
  GridKit::HiresComponent3<double, int>* comp3;

  Hires()
  {
    system   = new System();
    comp1    = new GridKit::HiresComponent1<double, int>(&bus, 1);
    bus_comp = new GridKit::HiresBus<double, int>(&bus, 2);
    comp3    = new GridKit::HiresComponent3<double, int>(&bus, 3);

    system->addComponent(comp1);
    system->addComponent(comp3);
    system->addComponent(bus_comp);
    system->addNode(&bus);
    system->allocate();

    bus.tag()[0] = true;
    bus.tag()[1] = true;
    system->tagDifferentiable();
    std::cout << "size tag " << bus_comp->tag()[0] << " " << bus_comp->tag()[1] << std::endl;

    for (auto i : system->tag())
    {
      std::cout << "tag " << i << std::endl;
    }

    // Initial state: y = (1, 0, 0, 0, 0, 0, 0, 0.0057) in HIRES numbering
    double* y  = system->y().getData();
    double* yp = system->yp().getData();
    for (int i = 0; i < system->size(); ++i)
    {
      y[i]  = 0.0;
      yp[i] = 0.0;
    }
    y[0] = 1.0;    // y1
    y[5] = 0.0057; // y8

    // Consistent y' = g(y): with y' = 0 the residual F = g(y) - y' equals g(y)
    system->updateTime(0.0, 0.0);
    system->evaluateResidual();
    const double* F = system->getResidual().getData();
    for (int i = 0; i < system->size(); ++i)
    {
      yp[i] = F[i];
    }

    system->y().setDataUpdated();
    system->yp().setDataUpdated();
  }

  ~Hires()
  {
    delete system;
  }
};

std::vector<double> solveMonolithic()
{
  Hires hires;

  Ida ida(hires.system);
  ida.setTolerance(1e-14, 1e-14);
  ida.setMaxSteps(1000000);
  ida.configureSimulation();
  ida.getDefaultInitialCondition();
  ida.initializeSimulation(0.0, false);
  ida.runSimulation(T_END);

  const double* y = hires.system->y().getData();
  return std::vector<double>(y, y + hires.system->size());
}

std::vector<double> solveCosim(int n_steps)
{
  Hires hires;

  // Partition 1: Component1 + Bus (+ interface carrying Component3's bus terms)
  // Partition 2: Component3
  auto* bus_interface = new Interface(&hires.bus, hires.comp3, 4);
  bus_interface->allocate();

  auto* partition1 = new Subsystem();
  partition1->addComponent(hires.comp1);
  partition1->addComponent(hires.bus_comp);
  partition1->addInterface(bus_interface);
  partition1->addNode(&hires.bus);

  auto* partition2 = new Subsystem();
  partition2->addComponent(hires.comp3);

  // Macro time points 0, dt, 2dt, ..., T_END
  std::vector<double> time_points(static_cast<size_t>(n_steps + 1));
  for (size_t k = 0; k <= static_cast<size_t>(n_steps); ++k)
  {
    time_points[k] = T_END * static_cast<double>(k) / static_cast<double>(n_steps);
  }

  size_t num_iterations = 1;
  size_t num_stages     = 4;

  auto* cosim = new Cosim(num_iterations, num_stages);
  cosim->addSubsystems(partition1);
  cosim->addSubsystems(partition2);
  cosim->allocate();
  cosim->initializeSimulation(time_points, &hires.system->y(), &hires.system->yp());
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

  const std::vector<int> step_counts = {200}; // {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5000};
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

  std::cout << "\n     dt     max error\n";
  for (size_t k = 0; k < step_counts.size(); ++k)
  {
    std::cout << std::scientific << std::setprecision(3)
              << std::setw(12) << T_END / step_counts[k] << ", "
              << std::setw(14) << errors[k]
              << std::endl;
  }

  return 0;
}
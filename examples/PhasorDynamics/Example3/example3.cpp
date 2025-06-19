#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>

// #include <sundials_core.h>
#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include "Example3_Powerworld_Reference.hpp"
#include "Model/PhasorDynamics/Branch/Branch.cpp"
#include "Model/PhasorDynamics/Branch/Branch.hpp"
#include "Model/PhasorDynamics/Bus/Bus.cpp"
#include "Model/PhasorDynamics/Bus/Bus.hpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.cpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.hpp"
#include "Model/PhasorDynamics/BusFault/BusFault.hpp"
#include "Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp"
#include "Model/PhasorDynamics/Governor/TGOV1/TurbineGov.hpp"
#include "Model/PhasorDynamics/SystemModel.hpp"
#include "Solver/Dynamic/Ida.cpp"
#include "Solver/Dynamic/Ida.hpp"
#include <Utilities/Testing.hpp>

#define _CRT_SECURE_NO_WARNINGS

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  printf("Example 3\n");

  /* Create model parts */
  SystemModel<double, size_t> sys;
  Bus<double, size_t>         bus1(0.9949877346411762, 0.09999703952427966);
  BusInfinite<double, size_t> bus2(1.0, 0.0);
  Branch<double, size_t>      branch(&bus1, &bus2, 0, 0.1, 0, 0);
  BusFault<double, size_t>    fault(&bus1, 0, 1e-3, 0);

  // Decleration
  Genrou<double, size_t>* gen;
  TurbineGov<double, size_t>* gov;

  // Instatiation
  gen = new Genrou<double, size_t>(&bus1,
                             1,
                             gov
                             1.,
                             0.05013,
                             3.,
                             0.,
                             0.,
                             7.,
                             .04,
                             .05,
                             .75,
                             2.1,
                             0.2,
                             0.18,
                             0.5,
                             0.5,
                             0.18,
                             0.15,
                             0.,
                             0.);

  // Governor of Generator
  gov = TurbineGov<double, size_t>(gen);


  /* Connect everything together */
  sys.addBus(&bus1);
  sys.addBus(&bus2);
  sys.addComponent(&branch);
  sys.addComponent(&fault);
  sys.addComponent(gen);
  sys.addComponent(gov);
  sys.allocate();

  double dt = 1.0 / 4.0 / 60.0;

  std::stringstream buffer;

  /* Set up simulation */
  Ida<double, size_t> ida(&sys);
  ida.configureSimulation();

  /* Run simulation */
  double start = static_cast<double>(clock());
  // ida.printOutputF(0, 0, buffer);
  ida.initializeSimulation(0.0, false);
  ida.runSimulationFixed(0.0, dt, 1.0, buffer);
  fault.setStatus(1);
  ida.initializeSimulation(1.0, false);
  ida.runSimulationFixed(1.0, dt, 1.1, buffer);
  fault.setStatus(0);
  ida.initializeSimulation(1.1, false);
  ida.runSimulationFixed(1.1, dt, 10.0, buffer);
  double stop = static_cast<double>(clock());

  // Go to the beginning of the data buffer
  buffer.seekg(0, std::ios::beg);

  double data;

  size_t i       = 0;   // data row counter
  size_t j       = 0;   // data column counter
  double Vr      = 0.0; // Bus real voltage
  double Vi      = 0.0; // Bus imaginary voltage
  double dw      = 0.0; // Generator frequency deviation [rad/s]
  double ti      = 0.0; // time
  double error_V = 0.0; // error in |V|

  // Read through the simulation data storred in the buffer
  while (buffer >> data)
  {
    // At the end of each data line compare computed data to Powerworld results
    // and reset column counter to zero.
    if ((i % 48) == 0)
    {
      double err =
          std::abs(std::sqrt(Vr * Vr + Vi * Vi) - reference_solution[i / 48][2])
          / (1.0 + std::abs(reference_solution[i / 48][2]));
      if (err > error_V)
        error_V = err;
      // std::cout << "t = " << ti << ": Vr = " << Vr << ", Vi = " << Vi << ", dw = " << dw;
      std::cout << "GridKit: t = " << ti
                << ", |V| = " << std::sqrt(Vr * Vr + Vi * Vi)
                << ", w = " << (1.0 + dw) << "\n";
      std::cout << "Ref    : t = " << reference_solution[i / 48][0]
                << ", |V| = " << reference_solution[i / 48][2]
                << ", w = " << reference_solution[i / 48][1]
                << "\n";
      std::cout << "Error in |V| = "
                << err
                << "\n";
      j  = 0;
      Vr = 0.0;
      Vi = 0.0;
      std::cout << "\n";
    }
    if (j == 0)
    {
      ti = data;
    }
    if (j == 2)
    {
      Vr = data;
    }
    if (j == 3)
    {
      Vi = data;
    }
    if (j == 5)
    {
      dw = data;
    }
    ++j;
    ++i;
    // if (i > 500)
    //   break;
  }

  // std::cout << buffer.str();
  int status = 0;
  std::cout << "Max error in |V| = " << error_V << "\n";
  if (error_V > 2e-4)
  {
    std::cout << "Test failed with error too large!\n";
    status = 1;
  }

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return 0;
}

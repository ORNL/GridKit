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

#include "Example1_Powerworld_Reference.hpp"
#include "Model/PhasorDynamics/Branch/Branch.cpp"
#include "Model/PhasorDynamics/Branch/Branch.hpp"
#include "Model/PhasorDynamics/Bus/Bus.cpp"
#include "Model/PhasorDynamics/Bus/Bus.hpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.cpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.hpp"
#include "Model/PhasorDynamics/BusFault/BusFault.hpp"
#include "Model/PhasorDynamics/SynchronousMachine/GENROUwS/GENROU.hpp"
#include "Model/PhasorDynamics/SystemModel.hpp"
#include "Solver/Dynamic/Ida.cpp"
#include "Solver/Dynamic/Ida.hpp"
#include <Utilities/Testing.hpp>

#define _CRT_SECURE_NO_WARNINGS

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  printf("Example 1 version 2\n");

  /* Create model parts */
  SystemModel<double, size_t> sys;
  Bus<double, size_t>         bus1(0.9949877346411762, 0.09999703952427966);
  BusInfinite<double, size_t> bus2(1.0, 0.0);
  Branch<double, size_t>      branch(&bus1, &bus2, 0, 0.1, 0, 0);
  BusFault                    fault(&bus1, 0, 1e-3, 0);

  GENROU gen(&bus1, 1, 1, 0.05013, 3, 0, 0, 7, .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0, 0);

  /* Connect everything together */
  sys.addBus(&bus1);
  sys.addBus(&bus2);
  sys.addComponent(&branch);
  sys.addComponent(&fault);
  sys.addComponent(&gen);
  sys.allocate();

  double dt = 1.0 / 4.0 / 60.0;

  /* Output file header */
  // FILE* f = fopen("example1_v2_results.csv", "w");
  // if (!f)
  //   printf("ERROR writing to output file!\n");
  // fprintf(f, "%s,%s", "t", "IDA Return Value");
  // for (int i = 0; i < sys.size(); ++i)
  //   fprintf(f, ",Y[%d]", i);
  // for (int i = 0; i < sys.size(); ++i)
  //   fprintf(f, ",Yp[%d]", i);
  // fprintf(f, "\n");

  std::stringstream buffer;
  // buffer << "t,IDA Return Value";
  // for (int i = 0; i < sys.size(); ++i)
  //   buffer << ",Y[" << i << "]";
  // for (int i = 0; i < sys.size(); ++i)
  //   buffer << ",Yp[" << i << "]";
  // buffer << std::endl;

  /* Set up simulation */
  Ida<double, size_t> ida(&sys);
  ida.configureSimulation();

  /* Run simulation */
  double start = (double) clock();
  // ida.printOutputF(0, 0, buffer);
  ida.initializeSimulation(0.0, false);
  ida.runSimulationFixed(0.0, dt, 1.0, buffer);
  fault.setStatus(1);
  ida.initializeSimulation(1.0, false);
  ida.runSimulationFixed(1.0, dt, 1.1, buffer);
  fault.setStatus(0);
  ida.initializeSimulation(1.1, false);
  ida.runSimulationFixed(1.1, dt, 10.0, buffer);

  buffer.seekg(0, std::ios::beg);
  double time;
  int    retvalue;
  double y0;

  // buffer >> time >> retvalue >> y0;
  // std::cout << time << " " << retvalue << " " << y0 << "\n";
  int    i  = 0;
  int    j  = 0;
  double Vr = 0.0;
  double Vi = 0.0;
  double dw = 0.0;
  double ti = 0.0;
  double error_V = 0.0;
  while (buffer >> time)
  {
    if ((i % 48) == 0)
    {
      double err = std::abs(std::sqrt(Vr * Vr + Vi * Vi) - reference_solution[i / 48][2]) / (1.0 + std::abs(reference_solution[i / 48][2]));
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
      ti = time;
    }
    else if (j == 1 + 1)
    {
      Vr = time;
    }
    else if (j == 2 + 1)
    {
      Vi = time;
    }
    else if (j == 4 + 1)
    {
      dw = time;
    }
    else
    {
      // do nothing
    }
    ++j;
    // std::cout << time << " ";
    ++i;
    // if (i > 500)
    //   break;
  }

  // std::cout << buffer.str();
  std::cout << "Max error in |V| = " << error_V << "\n";

  printf("\n\nComplete in %.4g seconds\n", (clock() - start) / CLOCKS_PER_SEC);
  // fclose(f);

  return 0;
}

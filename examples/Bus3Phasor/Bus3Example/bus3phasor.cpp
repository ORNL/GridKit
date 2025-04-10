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

#include "Model/PhasorDynamics/Branch/Branch.cpp"
#include "Model/PhasorDynamics/Branch/Branch.hpp"
#include "Model/PhasorDynamics/Bus/Bus.cpp"
#include "Model/PhasorDynamics/Bus/Bus.hpp"
#include "Model/PhasorDynamics/Load/Load.cpp"
#include "Model/PhasorDynamics/Load/Load.hpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.cpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.hpp"
#include "Model/PhasorDynamics/BusFault/BusFault.hpp"
#include "Model/PhasorDynamics/SynchronousMachine/GENROUwS/GENROU.hpp"
#include "Model/PhasorDynamics/SystemModel.hpp"
#include "Solver/Dynamic/Ida.cpp"
#include "Solver/Dynamic/Ida.hpp"

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
  Bus<double, size_t>         bus3(0.9949877346411762, 0.09999703952427966);
  Branch<double, size_t>      branch1(&bus1, &bus2, 0, 0.1, 0, 0);
  Branch<double, size_t>      branch2(&bus1, &bus3, 0, 0.1, 0, 0);
  Branch<double, size_t>      branch3(&bus3, &bus2, 0, 0.1, 0, 0);
  Load<double, size_t>        load1(&bus1, 0.2, 0.3);
  Load<double, size_t>        load2(&bus3, 0.2, 0.3);

  GENROU gen1(&bus1, 1, 1, 0.05013, 3, 0, 0, 7, .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0, 0);
  GENROU gen2(&bus1, 1, 1, 0.05013, 3, 0, 0, 7, .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0, 0);

  /* Connect everything together */
  sys.addBus(&bus1);
  sys.addBus(&bus2);
  sys.addBus(&bus3);
  sys.addComponent(&branch1);
  sys.addComponent(&branch2);
  sys.addComponent(&branch3);
  sys.addComponent(&load1);
  sys.addComponent(&load2);
  sys.addComponent(&gen1);
  sys.addComponent(&gen2);
  sys.allocate();

  double dt = 1.0 / 4.0 / 60.0;

  /* Output file header */
  FILE* f = fopen("example1_v2_results.csv", "w");
  if (!f)
    printf("ERROR writing to output file!\n");
  fprintf(f, "%s,%s", "t", "IDA Return Value");
  for (int i = 0; i < sys.size(); ++i)
    fprintf(f, ",Y[%d]", i);
  for (int i = 0; i < sys.size(); ++i)
    fprintf(f, ",Yp[%d]", i);
  fprintf(f, "\n");

  /* Set up simulation */
  Ida<double, size_t> ida(&sys);
  ida.configureSimulation();

  /* Run simulation */
  double start = (double) clock();
  ida.printOutputF(0, 0, f);
  ida.initializeSimulation(0.0, false);
  ida.runSimulationFixed(0.0, dt, 1.0, f);

  printf("Complete in %.4g seconds\n", (clock() - start) / CLOCKS_PER_SEC);
  fclose(f);

  return 0;
}

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
#include "Model/PhasorDynamics/SynchronousMachine/GENROUwS/GEN2.hpp"
#include "Model/PhasorDynamics/SystemModel.hpp"
#include "Solver/Dynamic/Ida.cpp"
#include "Solver/Dynamic/Ida.hpp"

#define _CRT_SECURE_NO_WARNINGS

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  printf("Example 1 version GENERATION 2\n");

  SystemModel<double, size_t> sys;
  Bus<double, size_t>         bus1(0.9949877346411762, 0.09999703952427966);
  BusInfinite<double, size_t> bus2(1.0, 0.0);
  Branch<double, size_t>      branch(&bus1, &bus2, 0.0, 0.1, 0, 0);
  GEN2 gen(&bus1, 1,3, 0, 0, 0.2, 0.1, 0.2);
  sys.allocate();


  /* Connect everything together */
  sys.addBus(&bus1);
  sys.addBus(&bus2);
  sys.addComponent(&branch);
  sys.addComponent(&gen);
  
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
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

#include "Example2_Powerworld_Reference.hpp"
#include "Model/PhasorDynamics/Branch/Branch.cpp"
#include "Model/PhasorDynamics/Branch/Branch.hpp"
#include "Model/PhasorDynamics/Bus/Bus.cpp"
#include "Model/PhasorDynamics/Bus/Bus.hpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.cpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.hpp"
#include "Model/PhasorDynamics/BusFault/BusFault.hpp"
#include "Model/PhasorDynamics/Load/Load.hpp"
#include "Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp"
#include "Model/PhasorDynamics/SystemModel.hpp"
#include "Solver/Dynamic/Ida.cpp"
#include "Solver/Dynamic/Ida.hpp"
#include <Utilities/Testing.hpp>

#define _CRT_SECURE_NO_WARNINGS

int main()
{
    using namespace GridKit::PhasorDynamics;
    using namespace AnalysisManager::Sundials;

    printf("Example 2 version 1\n");

    /* Create model parts */
    SystemModel<double, size_t> sys;
    BusInfinite<double, size_t> bus1(1.06, 0.0);
    Bus<double, size_t>         bus2(1.0599558398065716, -0.009675621941024773);
    Bus<double, size_t>         bus3(0.9610827543495831, -0.13122476630506485);
    Branch<double, size_t>      branch12(&bus1, &bus2, 0.05, 0.21, 0, 0.1);
    Branch<double, size_t>      branch13(&bus1, &bus3, 0.06, 0.15, 0, 0.12);
    Branch<double, size_t>      branch23(&bus2, &bus3, 0.08, 0.27, 0, 0.45);
    Genrou<double, size_t> gen2(&bus2, 1, 0.5, -0.07588, 2.7, 0., 0., 7., .04,
          .05, .75, 1.9, 0.17, 0.15, 0.4, 0.35, 0.15, 0.14999, 0., 0.);
    Genrou<double, size_t> gen3(&bus3, 1, 0.25, 0.26587, 1.6, 0., 0., 7.5, .04,
        .05, .75, 2.3, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
    //Genrou<double, size_t> gen3(&bus3, 1, 0.25, 0.26587, 1.6, 0, 0, 7.5, .04, 
    //    .05, .75, 2.3, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0, 0);
    //Genrou<double, size_t> gen3(&bus3, 1, 1.,  0.05013, 3., 0., 0., 7., .04,
    //    .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
    Load<double, size_t> load3(&bus3, 0.4447197839297772, 0.20330047265361242);
    BusFault<double, size_t>    fault(&bus3, 0, 1e-5, 0);

    /* Connect everything together */
    sys.addBus(&bus1);
    sys.addBus(&bus2);
    sys.addBus(&bus3);
    sys.addComponent(&branch12);
    sys.addComponent(&branch13);
    sys.addComponent(&branch23);
    sys.addComponent(&fault);
    sys.addComponent(&load3);
    sys.addComponent(&gen2);
    sys.addComponent(&gen3);
    sys.allocate();
    /*sys.initialize();
    sys.tagDifferentiable();
    sys.evaluateResidual();
    for (int i = 0; i< 46; ++i)
    {
        printf("%d %g %g %g\n", i, sys.y()[i], sys.yp()[i], sys.getResidual()[i]);
    }*/

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

  /* Check worst-case error */
  double worst_error = 0;
  double worst_error_time = 0;

  int i, j;
  const int stride = 94;
  const int nt = 2401;
  double results[stride];
  buffer.seekg(0, std::ios::beg);
  FILE *f = fopen("example2_results.csv", "w");
  fprintf(f, "Time,gen2speed,gen3speed,v2mag,v3mag\n");
  fprintf(f, "0,1,1,1,1\n");
  for (i = 0; i < nt-1; ++i)
  {
    for (j = 0; j < stride; ++j) buffer >> results[j];
    //for (j = 0; j < stride; ++j) printf("%d %d %g\n", i, j, results[j]);
    double t = results[0];
    double tref = reference_solution[i+1][0];
    //printf("Time GridKit %g PowerWorld %g\n", t, tref);
    double gen2speed = 1 + results[7];
    double gen2speed_ref = reference_solution[i+1][1];
    double gen3speed = 1 + results[28];
    double gen3speed_ref = reference_solution[i+1][2];
    double v2mag = sqrt(results[2]*results[2]+results[3]*results[3]);
    double v2mag_ref = reference_solution[i+1][4];
    double v3mag = sqrt(results[4]*results[4]+results[5]*results[5]);
    double v3mag_ref = reference_solution[i+1][5];
    fprintf(f, "%g,%g,%g,%g,%g\n", t, gen2speed, gen3speed, v2mag, v3mag);
    double err = fmax(fmax(fabs(gen2speed-gen2speed_ref), 
        fabs(gen3speed-gen3speed_ref)), fmax(fabs(v2mag-v2mag_ref), 
        fabs(v3mag-v3mag_ref)));
    if (err > worst_error)
    {
        worst_error = err;
        worst_error_time = t;
    }
    //printf("%g %g %g %g\n", gen2speed, gen3speed, v2mag, v3mag);
    //printf("%g %g %g %g\n", gen2speed_ref, gen3speed_ref, v2mag_ref, v3mag_ref);
  }
  fclose(f);

    printf("Worst error %g at time t=%g\n", worst_error, worst_error_time);

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

}
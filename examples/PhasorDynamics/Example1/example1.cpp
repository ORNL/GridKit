#include "example1.hpp"

#define _CRT_SECURE_NO_WARNINGS

using namespace GridKit::PhasorDynamics;
using namespace AnalysisManager::Sundials;

int main()
{
    printf("Example 1 version 2\n");

    /* Create model parts */
    SystemModel<double, size_t> sys;
    Bus<double, size_t> bus1(0.9949877346411762, 0.09999703952427966);
    BusInfinite<double, size_t> bus2(1.0, 0.0);
    Branch<double, size_t> branch(&bus1, &bus2, 0, 0.1, 0, 0);
    BusFault fault(&bus1, 0, 1e-3, 0);
    GENROU gen(&bus1, 1, 1, 0.05013, 3, 0, 0, 7, .04, .05, .75, 2.1, 0.2, 0.18, 
        0.5, 0.5, 0.18, 0.15, 0, 0);

    /* Connect everything together */
    sys.addBus(&bus1);
    sys.addBus(&bus2);
    sys.addComponent(&branch);
    sys.addComponent(&fault);
    sys.addComponent(&gen);
    sys.allocate();

    double dt = 1.0/4.0/60.0;

    /* Output file header */
    FILE *f = fopen("example1_v2_results.csv", "w");
    if (!f) printf("ERROR writing to output file!\n");
    fprintf(f, "%s,%s", "t", "IDA Return Value");
    for (int i = 0; i < sys.size(); ++i) fprintf(f, ",Y[%d]", i);
    for (int i = 0; i < sys.size(); ++i) fprintf(f, ",Yp[%d]", i);
    fprintf(f, "\n");

    /* Set up simulation */
    Ida<double, size_t> ida(&sys);
    ida.configureSimulation();

    /* Run simulation */
    double start = (double) clock();
    ida.printOutputF(0, 0, f);
    ida.initializeSimulation(0, false);
    ida.runSimulationFixed(0, dt, 1, f);
    fault.setStatus(1);
    ida.initializeSimulation(1, false);
    ida.runSimulationFixed(1, dt, 1.1, f);
    fault.setStatus(0);
    ida.initializeSimulation(1.1, false);
    ida.runSimulationFixed(1.1, dt, 30, f);

    printf("Complete in %.4g seconds\n", (clock() - start) / CLOCKS_PER_SEC);
    fclose(f);

    return 0;

}
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <time.h>

#include "Model/PhasorDynamics/Branch/Branch.hpp"
#include "Model/PhasorDynamics/Bus/Bus.hpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.hpp"
#include "Model/PhasorDynamics/BusFault/BusFault.hpp"
#include "Model/PhasorDynamics/Load/Load.hpp"
#include "Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassical.hpp"
#include "Model/PhasorDynamics/SystemModel.hpp"
#include "Solver/Dynamic/Ida.hpp"

#define _CRT_SECURE_NO_WARNINGS

int main(int argc, char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  std::cout << "Example 1 version GENERATION 2" << std::endl;

  // bus voltages
  double vr = 1.0;
  double vi = 0.0;

  // branch parameters
  double R = 0.0; // line series resistance
  double X = 0.1; // line series reactance
  double G = 0.0; // line shunt conductance
  double B = 0.0; // line shunt charging

  // Generator parameters
  double p0  = 1.0;     // real power output
  double q0  = 0.05013; // reactive power output
  double H   = 3.0;     // Initia constant
  double D   = 0.0;     // Damping coefficient
  double Ra  = 0.0;     // Winding resistance
  double Xdp = 0.2;     // Machine reactance parameter

  SystemModel<double, size_t>  sys;
  Bus<double, size_t>          bus1(vr, vi);
  BusInfinite<double, size_t>  bus2(vr, vi);
  Branch<double, size_t>       branch(&bus1, &bus2, R, X, G, B);
  GenClassical<double, size_t> gen(&bus1, 1, p0, q0, H, D, Ra, Xdp);

  /* Connect everything together */
  sys.addBus(&bus1);
  sys.addBus(&bus2);
  sys.addComponent(&branch);
  sys.addComponent(&gen);
  sys.allocate();
  sys.initialize();

  std::vector<std::vector<double>> outputData;

  // callback for outputting solution
  auto output_cb = [&](double t)
  {
    std::vector<double> yval;

    yval.push_back(t);
    for (auto val : sys.y())
    {
      yval.push_back(val);
    }

    outputData.push_back(yval);
  };

  output_cb(0);

  /* Set up simulation */
  Ida<double, size_t> ida(&sys);
  ida.configureSimulation();

  /* Run simulation */
  double start = static_cast<double>(clock());
  ida.initializeSimulation(0.0, false);
  size_t nout = 50;
  ida.runSimulation(1.0, nout, output_cb);

  // write solution to file if the user passes in a file name
  if (argc >= 2)
  {
    std::ofstream outfile(argv[1]);

    if (!outfile)
    {
      std::cout << "ERROR writing to output file!" << std::endl;
    }

    outfile << "t";
    for (int i = 0; i < sys.size(); ++i)
    {
      outfile << ",Y[" + std::to_string(i) + "]";
    }
    outfile << "\n";

    for (auto y_tstep : outputData)
    {
      outfile << y_tstep[0];
      for (int i = 1; i < y_tstep.size(); i++)
      {
        outfile << "," << y_tstep[i];
      }
      outfile << "\n";
    }

    outfile.close();
  }
  else
  {
    std::cout << "To print the solution, add a file name(eg:./gen2_example output.csv)" << std::endl;
  }

  printf("Complete in %.4g seconds\n", (clock() - start) / CLOCKS_PER_SEC);

  return 0;
}

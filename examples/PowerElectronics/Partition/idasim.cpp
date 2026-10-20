//

#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>

#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>
#include <examples/PowerElectronics/ExamplesHelper/SystemAssembler.hpp>

using index_type = size_t;
using real_type  = double;

struct GenDataIndex
{
  bool refframe;
  index_type  idx_y11;
  index_type  idx_y12;
  index_type  idx_y13;
  index_type  idx_y14;
  index_type  idx_y15;
};

std::vector<real_type> linspace(real_type start, real_type end, std::size_t num_points);
int runMonolithic(index_type num_intervals_1, index_type num_intervals_2, real_type tol, index_type N_size, std::string type);

constexpr bool WRITE_OUTPUT = true;

 // Those generators whose data we want to output
 const std::vector<size_t> generator_ids = {9, 10, 11, 99, 9999};

int main()
{

  std::string timing = "timing" ;

  if(WRITE_OUTPUT)
  { 
    timing = ""; //If output is to be written, we do not consider timing
  } 

  std::string info_type = "1e-11_ref_sol_" + timing;

  runMonolithic(5,  20, 1e-11, 5000, info_type);
  runMonolithic(10, 30, 1e-11, 5000, info_type);
  runMonolithic(15, 40, 1e-11, 5000, info_type);
  runMonolithic(20, 50, 1e-11, 5000, info_type);
  runMonolithic(25, 60, 1e-11, 5000, info_type);
  runMonolithic(30, 70, 1e-11, 5000, info_type);
  runMonolithic(35, 80, 1e-11, 5000, info_type);

  runMonolithic(0, 140, 1e-2, 5000, "1e-2_mono_" + timing);
  runMonolithic(0, 140, 1e-3, 5000, "1e-3_mono_" + timing);
  runMonolithic(0, 140, 1e-4, 5000, "1e-4_mono_" + timing);
  runMonolithic(0, 140, 1e-5, 5000, "1e-5_mono_" + timing);
  runMonolithic(0, 140, 1e-6, 5000, "1e-6_mono_" + timing);
  runMonolithic(0, 140, 1e-7, 5000, "1e-7_mono_" + timing);
  runMonolithic(0, 140, 1e-8, 5000, "1e-8_mono_" + timing);

  runMonolithic(0, 140, 1e-11, 5000, "1e-11_mono_reference" + timing);
}

int runMonolithic(index_type num_intervals_1, index_type num_intervals_2, real_type tol, index_type N_size, std::string type)
{
  using namespace GridKit;

  constexpr index_type SCALE_MICROGRID_MAX_STEPS = 100000;
  constexpr bool       use_jac                   = true;

  constexpr real_type t_init  = 0.0;
  constexpr real_type t_final = 1.0;
  constexpr size_t load_faulted_index = 10;

   real_type  rel_tol          = tol;
   real_type  abs_tol          = tol;


  // ---------------------------------------------------------------------------
  // Build the microgrid
  // ---------------------------------------------------------------------------

  ScaleMicrogridNetwork<real_type, index_type> network(N_size);

  auto* system = new PowerElectronicsModel<real_type, index_type>(use_jac);

  assembleSystem(network, *system);
  system->allocate();

  // ---------------------------------------------------------------------------
  // Collect the indices of the generators to inspect
  // ---------------------------------------------------------------------------
  std::vector<GenDataIndex> gen_data;

  for (const auto i : generator_ids)
  {
    const auto* gen = network.generators[i];
    const bool refframe = (i == 0);

    gen_data.push_back({
        refframe,
        gen->getNodeConnection(11),
        gen->getNodeConnection(12),
        gen->getNodeConnection(13),
        gen->getNodeConnection(14),
        refframe ? -1 : gen->getNodeConnection(15)
    });
  }

  // Apply a fault to the selected load
  for (size_t i = 0; i < network.loads.size(); i += 2)
  {
    if(i == load_faulted_index)
    {
      network.loads[i]->setResistance(3);
      network.loads[i]->setInductance(0.1);
    }
    else
    {
      network.loads[i]->setResistance(5);
      network.loads[i]->setInductance(20);
    }
  }

  // ---------------------------------------------------------------------------
  // Initial conditions
  // ---------------------------------------------------------------------------

  auto* y  = system->y().getData();
  auto* yp = system->yp().getData();

  for (index_type i = 0; i < system->size(); ++i)
  {
    y[i]  = 0.0;
    yp[i] = 0.0;
  }

  // Create initial derivatives specifics generated in MATLAB
  for (index_type i = 0; i < 2 * N_size; i++)
  {
    const auto& params = network.DGParam_list[i];

    const index_type offset = (i == 0) ? 0 : 13 * i - 1;

    yp[offset + 2] = params.Vn_;
    yp[offset + 4] = params.Kpv_ * params.Vn_;
    yp[offset + 6] = (params.Kpc_ * params.Kpv_ * params.Vn_) / params.Lf_;
  }

  // since the intial P_com = 0, the set the intial vector to the reference frame
  y[network.dg_signal.getNodeConnection(0).idx_] = network.DGParam_list[0].wb_;

  system->y().setDataUpdated();
  system->yp().setDataUpdated();

  // ---------------------------------------------------------------------------
  // Configure IDA
  // ---------------------------------------------------------------------------

  auto* idas = new AnalysisManager::Sundials::Ida<real_type, index_type>(system);

  idas->setTolerance(rel_tol, abs_tol);
  idas->setMaxSteps(SCALE_MICROGRID_MAX_STEPS);
  idas->configureSimulation();
  idas->getDefaultInitialCondition();
  idas->initializeSimulation(t_init);

  // ---------------------------------------------------------------------------
  // Run monolithic IDA
  // ---------------------------------------------------------------------------

  std::vector<real_type> Tspan1;
  std::vector<real_type> Tspan2;
  std::vector<real_type> Tspan;

    if (num_intervals_1 == 0)
  {
    Tspan = linspace(t_init, t_final, num_intervals_2 + 1);
  }
  else
  {
    Tspan1 = linspace(t_init, 0.01, num_intervals_1 + 1);
    Tspan2 = linspace(0.01, t_final, num_intervals_2 + 1);

    Tspan = Tspan1;
    Tspan.insert(Tspan.end(), Tspan2.begin() + 1, Tspan2.end());
  }


  const std::string filename = "ida_solution_intervals_"
                               + std::to_string(num_intervals_1)
                               + "_" + std::to_string(num_intervals_2)
                               + "_" + type 
                               + ".csv";

  std::ofstream output(filename);

  output << std::setprecision(16);

  if constexpr (WRITE_OUTPUT)
  {
    output << "time";

    for (const auto i : generator_ids)
    {
      output << ",delta_" << i
            << ",P_" << i
            << ",Q_" << i;
    }

    output << '\n';
  }

  auto out_cb = [&](const auto time) -> void
  {
    const auto* state_data = system->y().getData(GridKit::memory::HOST);

    output << time;

    for (const auto& gen : gen_data)
    {
      const auto y11 = state_data[gen.idx_y11];
      const auto y12 = state_data[gen.idx_y12];
      const auto y13 = state_data[gen.idx_y13];
      const auto y14 = state_data[gen.idx_y14];

      const auto delta = gen.refframe ? 0.0 : state_data[gen.idx_y15];
      const auto P     = y11 * y13 + y12 * y14;
      const auto Q     = -y11 * y14 + y12 * y13;

      output << "," << delta
            << "," << P
            << "," << Q;
    }

    output << '\n';
  };

  if constexpr (WRITE_OUTPUT)
  {
    out_cb(t_init);
  }

  const auto start = std::chrono::steady_clock::now();

  for (std::size_t i = 1; i < Tspan.size(); ++i)
  {
    idas->runSimulation(Tspan[i]);
    if constexpr (WRITE_OUTPUT)
    {
      out_cb(Tspan[i]);
    }
  }

  const auto                          end     = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = end - start;

  if constexpr (!WRITE_OUTPUT)
  {
    output << "\n# Metadata\n";
    output << "# solver: Monolithic IDA\n";
    output << "# rel_tol: " << rel_tol << '\n';
    output << "# abs_tol: " << abs_tol << '\n';
    output << "# num_intervals_1: " << num_intervals_1 << '\n';
    output << "# num_intervals_2: " << num_intervals_2 << '\n';
    output << "# num_output_steps: " << Tspan.size() - 1 << '\n';
    output << "# runtime_seconds: " << elapsed.count() << '\n';
  }
 
  output.close();

  std::cout << "\nIDA simulation time: "
            << elapsed.count() << " s\n";

  // ---------------------------------------------------------------------------
  // Clean up
  // ---------------------------------------------------------------------------

  delete idas;
  delete system;

  return 0;
}

std::vector<real_type> linspace(real_type start, real_type end, std::size_t num_points)
{
  if (num_points == 0)
  {
    return {};
  }

  if (num_points == 1)
  {
    return {start};
  }

  std::vector<real_type> points(num_points);

  const real_type step = (end - start) / static_cast<real_type>(num_points - 1);

  for (std::size_t i = 0; i < num_points; ++i)
  {
    points[i] = start + static_cast<real_type>(i) * step;
  }

  // Guarantee the endpoint despite floating-point roundoff.
  points.back() = end;

  return points;
}

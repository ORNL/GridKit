// 

#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>
#include <examples/PowerElectronics/ExamplesHelper/SystemAssembler.hpp>

using index_type = size_t;
using real_type  = double;

std::vector<real_type> linspace(real_type start, real_type end, std::size_t num_points);


struct GenDataIndex
{
  bool       refframe;
  index_type idx_y8;
  index_type idx_y9;
  index_type idx_y10;
  index_type idx_y11;
  index_type idx_y12;
};

int main()
{
  using namespace GridKit;

  constexpr index_type SCALE_MICROGRID_MAX_STEPS = 1000000;
  constexpr index_type N_size  = 5000;
  constexpr bool       use_jac = true;

  constexpr real_type t_init  = 0.0;
  constexpr real_type t_final = 1.0;

  constexpr real_type  rel_tol         = 1e-10;
  constexpr real_type  abs_tol         = 1e-10;
  constexpr index_type num_intervals_1 = 35;
  constexpr index_type num_intervals_2 = 70;

  // ---------------------------------------------------------------------------
  // Build the microgrid
  // ---------------------------------------------------------------------------

  ScaleMicrogridNetwork<real_type, index_type> network(N_size);

  auto* system = new PowerElectronicsModel<real_type, index_type>(use_jac);

  assembleSystem(network, *system);
  system->allocate();

  // ---------------------------------------------------------------------------
  // Identify generator variables used to compute delta, P, and Q
  // ---------------------------------------------------------------------------

  std::vector<GenDataIndex> genData(5);

  for (index_type i = 0; i < genData.size(); ++i)
  {
    auto* gen     = network.generators[i];
    auto& gendata = genData[i];

    gendata.refframe = (i == 0);
    gendata.idx_y8   = gen->getNodeConnection(8);
    gendata.idx_y9   = gen->getNodeConnection(9);
    gendata.idx_y10  = gen->getNodeConnection(10);
    gendata.idx_y11  = gen->getNodeConnection(11);
    gendata.idx_y12  = gendata.refframe ? 0 : gen->getNodeConnection(12);
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

  for (index_type i = 0; i < 2 * N_size; ++i)
  {
    const auto& params = network.DGParam_list[i];

    yp[13 * i - 1 + 2] = params.Vn_;
    yp[13 * i - 1 + 4] = params.Kpv_ * params.Vn_;
    yp[13 * i - 1 + 6] = (params.Kpc_ * params.Kpv_ * params.Vn_) / params.Lf_;
  }

  // Initial P_com = 0, so initialize the reference-frame frequency.
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

  std::vector<real_type> time_points = linspace(0.0, 0.01, num_intervals_1 + 1);

  std::vector<real_type> time_points1 = linspace(0.01, 1.0, num_intervals_2 + 1);

  time_points.insert(time_points.end(), time_points1.begin() + 1, time_points1.end());

  const std::string filename = "ida_solution_intervals_" + std::to_string(num_intervals_1) +
                               "_" + std::to_string(num_intervals_2) + "_ref.csv";

  std::ofstream output(filename);

  output << std::setprecision(10);

  auto out_cb = [&](const auto time) -> void
  {
    const auto* state_data = system->y().getData(GridKit::memory::HOST);

    output << std::scientific << std::setprecision(8)
           << std::setw(18) << time;

    for (index_type i = 0; i < genData.size(); ++i)
    {
      const auto& gendata = genData[i];

      const auto y8  = state_data[gendata.idx_y8];
      const auto y9  = state_data[gendata.idx_y9];
      const auto y10 = state_data[gendata.idx_y10];
      const auto y11 = state_data[gendata.idx_y11];

      const auto delta = gendata.refframe ? 0.0 : state_data[gendata.idx_y12];

      const auto P = y8 * y10 + y9 * y11;
      const auto Q = -y8 * y11 + y9 * y10;

      output << "; "
             << std::setw(12) << delta << ", "
             << std::setw(12) << P     << ", "
             << std::setw(12) << Q;
    }

    output << '\n';
  };

  out_cb(t_init);

  const auto start = std::chrono::steady_clock::now();

  for (std::size_t i = 1; i < time_points.size(); ++i)
  {
    idas->runSimulation(time_points[i]);
    out_cb(time_points[i]);
  }

  const auto end = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = end - start;

  output << "\n# Metadata\n";
  output << "# solver: Monolithic IDA\n";
  output << "# rel_tol: " << rel_tol << '\n';
  output << "# abs_tol: " << abs_tol << '\n';
  output << "# num_intervals_1: " << num_intervals_1 << '\n';
  output << "# num_intervals_2: " << num_intervals_2 << '\n';
  output << "# num_output_steps: " << time_points.size() - 1 << '\n';
  output << "# runtime_seconds: " << elapsed.count() << '\n';

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

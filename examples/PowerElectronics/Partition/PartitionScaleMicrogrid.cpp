#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/SignalNode.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

using index_type = size_t;
using real_type  = double;

int printMicrogridSystems(index_type N_size);

/**
 * @brief Run Scale Microgrid for a specific N given by the user.
 *
 * @param[in] argc should be 1 if no N given, 2 if N given
 * @param[in] argv can contain the N value (argv[1])
 * @return int
 */
int main(int argc, char const* argv[])
{
  // Default value
  index_type N_size = 2;

  // Parse command line arguments if provided
  if (argc > 1)
  {
    try
    {
      N_size = static_cast<index_type>(std::stoi(argv[1]));
    }
    catch (const std::exception& e)
    {
      std::cerr << "Error parsing grid size argument: " << e.what() << std::endl;
      std::cerr << "Using default grid size = " << N_size << std::endl;
    }
  }

  return printMicrogridSystems(N_size);
}

/**
 * @brief Tests network of distributed generators.
 *
 * @param[in] N_size - The number of DG line load cobinations to generate for scale
 * @return int returns 0 if successful, >0 otherwise
 */
int printMicrogridSystems(index_type N_size)
{
  using namespace GridKit;

  const size_t num_ibrs       = 2 * N_size;
  const size_t num_partitions = 4;

  assert(num_partitions <= num_ibrs);

  bool use_jac = true;

  real_type rel_tol = 1e-5;
  real_type abs_tol = 1e-5;

  // Create circuit model
  auto* sys_model = new PowerElectronicsModel<real_type, index_type>(rel_tol,
                                                                     abs_tol,
                                                                     use_jac,
                                                                     2000);

  // Ensure minimum size requirement
  if (N_size < 1)
  {
    std::cout << "N_size must be at least 1.\n";
    return 1;
  }

  // Modeled after the problem in the paper
  // Every Bus has the same virtual resistance. This is due to numerical stability as mentioned in the paper.
  real_type RN = 1.0e4;

  // DG Params Vector
  // All DGs have the same set of parameters except for the first two.
  GridKit::DistributedGeneratorParameters<real_type, index_type> DG_parms1;
  DG_parms1.wb_  = 2.0 * M_PI * 50.0;
  DG_parms1.wc_  = 31.41;
  DG_parms1.mp_  = 9.4e-5;
  DG_parms1.Vn_  = 380.0;
  DG_parms1.nq_  = 1.3e-3;
  DG_parms1.F_   = 0.75;
  DG_parms1.Kiv_ = 420.0;
  DG_parms1.Kpv_ = 0.1;
  DG_parms1.Kic_ = 2.0e4;
  DG_parms1.Kpc_ = 15.0;
  DG_parms1.Cf_  = 5.0e-5;
  DG_parms1.rLf_ = 0.1;
  DG_parms1.Lf_  = 1.35e-3;
  DG_parms1.rLc_ = 0.03;
  DG_parms1.Lc_  = 0.35e-3;

  GridKit::DistributedGeneratorParameters<real_type, index_type> DG_parms2;
  DG_parms2.wb_  = 2.0 * M_PI * 50.0;
  DG_parms2.wc_  = 31.41;
  DG_parms2.mp_  = 12.5e-5;
  DG_parms2.Vn_  = 380.0;
  DG_parms2.nq_  = 1.5e-3;
  DG_parms2.F_   = 0.75;
  DG_parms2.Kiv_ = 390.0;
  DG_parms2.Kpv_ = 0.05;
  DG_parms2.Kic_ = 16.0e3;
  DG_parms2.Kpc_ = 10.5;
  DG_parms2.Cf_  = 50.0e-6;
  DG_parms2.rLf_ = 0.1;
  DG_parms2.Lf_  = 1.35e-3;
  DG_parms2.rLc_ = 0.03;
  DG_parms2.Lc_  = 0.35e-3;

  std::vector<GridKit::DistributedGeneratorParameters<real_type, index_type>> DGParams_list(2 * N_size, DG_parms2);

  // First two generators use parameters 1
  if (DGParams_list.size() >= 1)
    DGParams_list[0] = DG_parms1;
  if (DGParams_list.size() >= 2)
    DGParams_list[1] = DG_parms1;

  // line vector params
  // Every odd line has the same parameters and every even line has the same parameters
  real_type              rline1 = 0.23;
  real_type              Lline1 = 0.1 / (2.0 * M_PI * 50.0);
  real_type              rline2 = 0.35;
  real_type              Lline2 = 0.58 / (2.0 * M_PI * 50.0);
  std::vector<real_type> rline_list(2 * N_size - 1, 0.0);
  std::vector<real_type> Lline_list(2 * N_size - 1, 0.0);
  for (index_type i = 0; i < rline_list.size(); i++)
  {
    rline_list[i] = (i % 2) ? rline2 : rline1;
    Lline_list[i] = (i % 2) ? Lline2 : Lline1;
  }

  // load parms
  // Only the first load has the same paramaters.
  real_type rload1 = 3.0;
  real_type Lload1 = 2.0 / (2.0 * M_PI * 50.0);
  real_type rload2 = 2.0;
  real_type Lload2 = 1.0 / (2.0 * M_PI * 50.0);

  std::vector<real_type> rload_list(N_size, rload2);
  std::vector<real_type> Lload_list(N_size, Lload2);
  if (rload_list.size() >= 1)
  {
    rload_list[0] = rload1;
    Lload_list[0] = Lload1;
  }

  using SignalNode = GridKit::PowerElectronics::SignalNode<real_type, index_type>;
  using Bus        = GridKit::PowerElectronics::MicrogridBus<real_type, index_type>;
  using BusDQ      = MicrogridBusDQ<real_type, index_type>;
  using Generator  = DistributedGenerator<real_type, index_type>;
  using Line       = MicrogridLine<real_type, index_type>;
  using Load       = MicrogridLoad<real_type, index_type>;

  std::vector<Bus*>       buses(num_ibrs, nullptr);
  std::vector<BusDQ*>     busesDQ(num_ibrs, nullptr);
  std::vector<Generator*> generators(num_ibrs, nullptr);
  std::vector<Line*>      lines(num_ibrs, nullptr);
  std::vector<Load*>      loads(num_ibrs, nullptr);

  SignalNode dg_signal;
  sys_model->addNode(&dg_signal);

  for (size_t i = 0; i < 2 * N_size; i++)
  {
    buses[i] = new Bus();
    sys_model->addNode(buses[i]);
  }

  // Create the reference DG
  auto* dg_ref = new DistributedGenerator<real_type, index_type>(0,
                                                                 DGParams_list[0],
                                                                 true,
                                                                 &dg_signal,
                                                                 buses[0]);
  sys_model->addComponent(dg_ref);

  generators[0] = dg_ref;

  // Keep track of models and index location
  index_type model_id = 1;
  // Add all other DGs
  for (index_type i = 1; i < 2 * N_size; i++)
  {
    // current DG to add
    auto* dg = new DistributedGenerator<real_type, index_type>(model_id++,
                                                               DGParams_list[i],
                                                               false,
                                                               &dg_signal,
                                                               buses[i]);

    generators[i] = dg;
    sys_model->addComponent(dg);
  }

  // Load all the Line components
  for (index_type i = 0; i < 2 * N_size - 1; i++)
  {
    // line
    auto* line_model = new MicrogridLine<real_type, index_type>(model_id++,
                                                                rline_list[i],
                                                                Lline_list[i],
                                                                &dg_signal,
                                                                buses[i],
                                                                buses[i + 1]);
    lines[i + 1]     = line_model;
    sys_model->addComponent(line_model);
  }

  //  Load all the Load components
  for (index_type i = 0; i < N_size; i++)
  {
    auto* load_model = new MicrogridLoad<real_type, index_type>(model_id++,
                                                                rload_list[i],
                                                                Lload_list[i],
                                                                &dg_signal,
                                                                buses[2 * i]);
    loads[2 * i]     = load_model;
    sys_model->addComponent(load_model);
  }

  // Add all the microgrid Virtual DQ Buses
  for (index_type i = 0; i < 2 * N_size; i++)
  {
    auto* virDQbus_model = new MicrogridBusDQ<real_type, index_type>(model_id++, RN, buses[i]);

    busesDQ[i] = virDQbus_model;
    sys_model->addComponent(virDQbus_model);
  }

  sys_model->allocate();

  std::vector<double> y;
  std::vector<double> yp;

  for (size_t i = 0; i < sys_model->size(); i++)
  {
    y.push_back(static_cast<double>(i + 1));
    yp.push_back(static_cast<double>(i + 1));
  }

  for (size_t i = 0; i < sys_model->size(); i++)
  {
    sys_model->y()[i]  = y[i];
    sys_model->yp()[i] = yp[i];
  }

  sys_model->evaluateResidual();
  std::vector<double> f_sysmodel = sys_model->getResidual();

  std::vector<SubsystemModel<real_type, index_type>*> subsystems(num_partitions);

  // Partition the system
  index_type q     = (num_ibrs) / num_partitions;
  index_type r     = (num_ibrs) % num_partitions;
  index_type index = 0;
  for (index_type j = 0; j < num_partitions; j++)
  {
    auto* partition = new SubsystemModel<real_type, index_type>(false);

    // add Reference rotor to the first partition
    if (j == 0)
    {
      partition->addNode(&dg_signal);
    }

    std::cout << "Partition " << j << std::endl;

    index_type part_size = q + (j < r ? 1 : 0);
    index_type end       = std::min(index + part_size, num_ibrs);

    // Add all components belonging to this partition
    for (; index < end; ++index)
    {
      partition->addComponent(generators[index]);
      std::cout << "Comp Gen" << index << std::endl;
      partition->addComponent(busesDQ[index]);
      std::cout << "Comp BDQ" << index << std::endl;

      if (loads[index] != nullptr)
      {
        partition->addComponent(loads[index]);
        std::cout << "Comp Load" << index << std::endl;
      }

      if (lines[index] != nullptr)
      {
        partition->addComponent(lines[index]);
        std::cout << "Comp Line" << index << std::endl;
      }

      std::cout << "Comp Bus " << index << std::endl;

      partition->addNode(buses[index]);
    }

    // Add partition interface at a partition point
    if (index < num_ibrs)
    {
      auto* linecopy = new GridKit::MicrogridLine<real_type, index_type>(*lines[index]);

      auto* busInterface = new GridKit::BusPartitionInterface<double, size_t>(*buses[index - 1], *linecopy, model_id++);
      busInterface->allocate();

      std::cout << "Interface " << "(" << index - 1 << ", " << index << ")" << std::endl;

      partition->addComponent(busInterface);
    }

    subsystems[j] = partition;
  }

  for (auto* partition : subsystems)
  {
    partition->toString();
    partition->allocate();
  }

  // Distribute externals to partition 1
  for (auto* partition : subsystems)
  {
    for (size_t i = 0; i < partition->getExternSize(); i++)
    {
      partition->getExternalDataY()[i]  = y[partition->getExternalIndices()[i]];
      partition->getExternalDataYP()[i] = yp[partition->getExternalIndices()[i]];
    }
  }

  for (auto* partition : subsystems)
  {
    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      partition->y()[i]  = y[partition->getNodeConnection(i)];
      partition->yp()[i] = yp[partition->getNodeConnection(i)];
    }
  }

  for (auto* partition : subsystems)
  {
    partition->evaluateResidual();
  }

  std::vector<double> f(sys_model->size(), 0.0);
  std::vector<double> error(sys_model->size(), 1.0);

  // Get internal residuals from partition 1
  for (auto* partition : subsystems)
  {
    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      f[partition->getNodeConnection(i)] = partition->getResidual()[i];
    }
  }

  for (size_t i = 0; i < sys_model->size(); i++)
  {
    error[i] = f_sysmodel[i] - f[i];
    std::cout << i << " " << error[i] << " ---------- " << i << std::endl;
  }

  double max_error = 0;
  for (size_t i = 0; i < sys_model->size(); i++)
  {
    if (max_error < std::abs(error[i]))
    {
      max_error = std::abs(error[i]);
    }
  }

  std::cout << "\nMax Error of Reference and Partition Evaluation: " << max_error << std::endl;

  return 0;
}

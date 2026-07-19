#include <omp.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>

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

#include "jac_test_helper.hpp"

/******************************************************************************
 * Partitioned Residual Evaluation
 *
 * Construct subsystem models by partitioning the original microgrid into
 * independent groups of components. Each subsystem is allocated and receives
 * the appropriate subset of the global state vectors (y and yp), consisting
 * of:
 *
 *   - External variables (coupling variables from neighboring partitions)
 *   - Internal variables (states owned by the partition)
 *
 * After distributing the state information, each partition independently
 * evaluates its residual. The local residuals are then scattered back into
 * the global residual vector to reconstruct the monolithic residual.
 *
 * Finally, the reconstructed residual is compared against the reference
 * monolithic evaluation to verify the correctness of the partitioning
 * implementation.
 ******************************************************************************/

using index_type = size_t;
using real_type  = double;

struct RunResult
{
  bool        success = true;
  index_type  num_partitions;
  real_type   partition_time;  // seconds
  real_type   monolithic_time; // seconds
  real_type   speedup;         // monolithic / partition
  real_type   max_error;
  std::string subsys_jac;
};

RunResult printMicrogridSystems(index_type N_size, index_type num_partitions);

/**
 * @brief Run Scale Microgrid for a specific N given by the user.
 *
 * @param[in] argc should be 1 if no N given, 2 if N given
 * @param[in] argv can contain the N value (argv[1])
 * @return int
 */
int main()
{
  index_type N_size = 5000;

  std::cout << std::format("{:<16}{:>16}{:>18}{:>12}{:>14}{:>16}\n",
                           "num_partitions",
                           "partition_time",
                           "monolithic_time",
                           "speedup",
                           "error",
                           "Jacobians");

  std::cout << std::string(93, '-') << "\n";

  for (index_type p : {10, 48, 100, 500, 1000, 2000, 5000})
  {
    assert(p <= 2 * N_size);

    RunResult r = printMicrogridSystems(N_size, p);

    if (!r.success)
    {
      return 1;
    }

    std::cout << std::format("{:<16d}{:>14.4f} s{:>16.4f} s{:>11.2f}x{:>14.3e} {:>16s}\n",
                             r.num_partitions,
                             r.partition_time,
                             r.monolithic_time,
                             r.speedup,
                             r.max_error,
                             r.subsys_jac);
  }

  return 0;
}

/**
 * @brief Tests network of distributed generators.
 *
 * @param[in] N_size - The number of DG line load cobinations to generate for scale
 * @return int returns 0 if successful, >0 otherwise
 */
RunResult printMicrogridSystems(index_type N_size, index_type num_partitions)
{
  using namespace GridKit;

  const index_type num_ibrs = 2 * N_size;

  assert(num_partitions <= num_ibrs);

  bool use_jac = true;

  // Create circuit model
  auto* sys_model_control = new PowerElectronicsModel<real_type, index_type>(use_jac);

  // Ensure minimum size requirement
  assert(N_size >= 1);

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
  {
    DGParams_list[0] = DG_parms1;
  }
  if (DGParams_list.size() >= 2)
  {
    DGParams_list[1] = DG_parms1;
  }

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

  using SignalNode         = GridKit::PowerElectronics::SignalNode<real_type, index_type>;
  using Bus                = GridKit::PowerElectronics::MicrogridBus<real_type, index_type>;
  using BusDQ              = MicrogridBusDQ<real_type, index_type>;
  using Generator          = DistributedGenerator<real_type, index_type>;
  using Line               = MicrogridLine<real_type, index_type>;
  using Load               = MicrogridLoad<real_type, index_type>;
  using PartitionInterface = BusPartitionInterface<real_type, index_type>;

  std::vector<Bus*>                buses(num_ibrs, nullptr);
  std::vector<BusDQ*>              busesDQ(num_ibrs, nullptr);
  std::vector<Generator*>          generators(num_ibrs, nullptr);
  std::vector<Line*>               lines(num_ibrs, nullptr);
  std::vector<Load*>               loads(num_ibrs, nullptr);
  std::vector<PartitionInterface*> partitionInterface;
  std::vector<Line*>               linesCopies;

  SignalNode dg_signal;
  // sys_model_control->addNode(&dg_signal);

  for (size_t i = 0; i < 2 * N_size; i++)
  {
    buses[i] = new Bus();
    // sys_model_control->addNode(buses[i]);
  }

  // Create the reference DG
  auto* dg_ref = new DistributedGenerator<real_type, index_type>(0,
                                                                 DGParams_list[0],
                                                                 true,
                                                                 &dg_signal,
                                                                 buses[0]);
  // sys_model_control->addComponent(dg_ref);

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
    // sys_model_control->addComponent(dg);
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
    // sys_model_control->addComponent(line_model);
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
    // sys_model_control->addComponent(load_model);
  }

  // Add all the microgrid Virtual DQ Buses
  for (index_type i = 0; i < 2 * N_size; i++)
  {
    auto* virDQbus_model = new MicrogridBusDQ<real_type, index_type>(model_id++, RN, buses[i]);

    busesDQ[i] = virDQbus_model;
    // sys_model_control->addComponent(virDQbus_model);
  }

  sys_model_control->addNode(&dg_signal);
  for (index_type i = 0; i < num_ibrs; ++i)
  {
    sys_model_control->addComponent(generators[i]);
    sys_model_control->addComponent(busesDQ[i]);

    if (loads[i] != nullptr)
    {
      sys_model_control->addComponent(loads[i]);
    }

    if (lines[i] != nullptr)
    {
      sys_model_control->addComponent(lines[i]);
    }
    sys_model_control->addNode(buses[i]);
  }

  sys_model_control->allocate();

  std::vector<real_type> y;
  std::vector<real_type> yp;

  for (size_t i = 0; i < sys_model_control->size(); i++)
  {
    y.push_back(static_cast<real_type>(i + 1));
    yp.push_back(static_cast<real_type>(i + 1));
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  auto* sys_model_y  = sys_model_control->y().getData();
  auto* sys_model_yp = sys_model_control->yp().getData();

  for (size_t i = 0; i < sys_model_control->size(); i++)
  {
    sys_model_y[i]  = y[i];
    sys_model_yp[i] = yp[i];
  }

  sys_model_control->y().setDataUpdated();
  sys_model_control->yp().setDataUpdated();

  sys_model_control->updateTime(2, 5);
  sys_model_control->evaluateResidual();

  auto end_time        = std::chrono::high_resolution_clock::now();
  auto monolithic_time = std::chrono::duration<real_type>(end_time - start_time);

  sys_model_control->evaluateJacobian();

  auto* f_sysmodel_control = sys_model_control->getResidual().getData();
  auto* full_jac           = sys_model_control->getCsrJacobian();

  //---------------------------------------------------------------
  // Partition the monolithic network into independent subsystem
  // models. Each subsystem owns a contiguous block of generators,
  // buses, lines, and loads. Whenever a partition boundary is
  // encountered, a BusPartitionInterface is inserted to expose the
  // neighboring partition's boundary variables while maintaining
  // the same electrical connectivity as the original network.
  //---------------------------------------------------------------

  std::vector<SubsystemModel<real_type, index_type>*> subsystems(num_partitions);

  // Partition the system
  index_type q     = (num_ibrs) / num_partitions;
  index_type r     = (num_ibrs) % num_partitions;
  index_type index = 0;
  for (index_type j = 0; j < num_partitions; j++)
  {
    auto* partition = new SubsystemModel<real_type, index_type>();

    // add Reference rotor to the first partition
    if (j == 0)
    {
      partition->addNode(&dg_signal);
    }

    index_type part_size = q + (j < r ? 1 : 0);
    index_type end       = std::min(index + part_size, num_ibrs);

    // Add all components belonging to this partition
    for (; index < end; ++index)
    {
      partition->addComponent(generators[index]);
      partition->addComponent(busesDQ[index]);

      if (loads[index] != nullptr)
      {
        partition->addComponent(loads[index]);
      }

      if (lines[index] != nullptr)
      {
        partition->addComponent(lines[index]);
      }

      partition->addNode(buses[index]);
    }
    // Add the partition interface to the left partition
    if (index < num_ibrs)
    {
      auto* linecopy     = new GridKit::MicrogridLine<real_type, index_type>(*lines[index]);
      auto* busInterface = new GridKit::BusPartitionInterface<real_type, index_type>(*buses[index - 1], *linecopy, model_id++);

      busInterface->allocate();
      partitionInterface.push_back(busInterface);
      linesCopies.push_back(linecopy);

      partition->addComponent(busInterface);
    }

    subsystems[j] = partition;
  }

  //---------------------------------------------------------------
  // Allocate each subsystem and evaluate the partitioned residual.
  //
  // The global state vectors (y and yp) are scattered into the
  // corresponding internal and external vectors of every subsystem.
  // After evaluating the residuals independently, the local
  // residuals are gathered back into the global residual vector for
  // comparison with the monolithic evaluation.
  //---------------------------------------------------------------

  std::vector<real_type> f(sys_model_control->size(), 1.0);
  std::vector<real_type> error(sys_model_control->size(), 1.0);

  for (auto* partition : subsystems)
  {
    partition->allocate();
    partition->updateTime(2, 5);
  }

  start_time = std::chrono::high_resolution_clock::now();

// Distribute externals to partition 1
#pragma omp parallel for schedule(guided)
  for (auto* partition : subsystems)
  {
    // Distribute external variables
    for (size_t i = 0; i < partition->getExternSize(); i++)
    {
      partition->getExternalDataY()[i]  = y[partition->getExternalIndices()[i]];
      partition->getExternalDataYP()[i] = yp[partition->getExternalIndices()[i]];
    }

    auto* partition_y  = partition->y().getData();
    auto* partition_yp = partition->yp().getData();

    // Distribute internal variables
    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      partition_y[i]  = y[partition->getNodeConnection(i)];
      partition_yp[i] = yp[partition->getNodeConnection(i)];
    }

    partition->y().setDataUpdated();
    partition->yp().setDataUpdated();

    // Evaluate residual of this partition
    partition->evaluateResidual();

    auto* residual = partition->getResidual().getData();
    // Reconstructs the monolithic residual from the partition residuals
    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      f[partition->getNodeConnection(i)] = residual[i];
    }
  }

  end_time            = std::chrono::high_resolution_clock::now();
  auto partition_time = std::chrono::duration<real_type>(end_time - start_time);

  //---------------------------------------------------------------
  // Compare the reconstructed partition residual against the
  // reference monolithic residual. The per-entry error is computed
  // and the maximum error is reported to verify that the partitioned
  // formulation reproduces the monolithic model within numerical
  // roundoff.
  //---------------------------------------------------------------

  bool matched = true;

  for (auto* partition : subsystems)
  {
    partition->evaluateJacobian();

    matched &= GridKit::Testing::verifySubsystemJacobian(
        *full_jac,
        *partition->getCsrJacobian(),
        *partition);
  }

  real_type max_error = 0;
  for (size_t i = 0; i < sys_model_control->size(); i++)
  {
    error[i] = abs(f[i] - f_sysmodel_control[i]) / (f_sysmodel_control[i] + 1);

    if (max_error < error[i])
    {
      max_error = error[i];
    }
  }

  RunResult result;
  result.num_partitions  = num_partitions;
  result.partition_time  = partition_time.count();
  result.monolithic_time = monolithic_time.count();
  result.speedup         = monolithic_time.count() / partition_time.count();
  result.max_error       = max_error;
  result.subsys_jac      = matched ? "Correct" : "Wrong";

  if (!matched)
  {
    std::cout << "ERROR: At least one subsystem Jacobian is incorrect!" << std::endl;
    result.success = false;
  }

  if (max_error > std::numeric_limits<double>::epsilon())
  {
    std::cout << "ERROR: Max Error too high!" << std::endl;
    result.success = false;
  }

  for (auto* linescpy : linesCopies)
  {
    delete linescpy;
  }
  for (auto* partition : subsystems)
  {
    delete partition;
  }
  for (auto* part_interface : partitionInterface)
  {
    delete part_interface;
  }
  delete sys_model_control;

  return result;
}

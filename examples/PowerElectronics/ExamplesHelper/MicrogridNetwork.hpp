// MicrogridNetwork.hpp

#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/SignalNode.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>

/*
 * Contains components and nodes that make up the scaled microgrid network.
 *
 * The network contains 2 * N_size IBRs and stores the physical components
 * that make up the scale microgrid.
 */
template <class ScalarT, typename IdxT>
struct ScaleMicrogridNetwork
{

  using SignalNode  = GridKit::PowerElectronics::SignalNode<ScalarT, IdxT>;
  using Bus         = GridKit::PowerElectronics::MicrogridBus<ScalarT, IdxT>;
  using BusDQ       = GridKit::MicrogridBusDQ<ScalarT, IdxT>;
  using DGGenerator = GridKit::DistributedGenerator<ScalarT, IdxT>;
  using Line        = GridKit::MicrogridLine<ScalarT, IdxT>;
  using Load        = GridKit::MicrogridLoad<ScalarT, IdxT>;
  using GenParams   = GridKit::DistributedGeneratorParameters<ScalarT, IdxT>;

  size_t                    model_id_next;
  size_t                    N_size;
  SignalNode                dg_signal;
  std::vector<Bus>          buses;
  std::vector<BusDQ*>       busesDQ;
  std::vector<DGGenerator*> generators;
  std::vector<Line*>        lines;
  std::vector<Load*>        loads;
  std::vector<GenParams>    DGParam_list;

  ScaleMicrogridNetwork(size_t n_size)
    : model_id_next(0),
      N_size(n_size),
      buses(2 * n_size),
      busesDQ(2 * n_size, nullptr),
      generators(2 * n_size, nullptr),
      lines(2 * n_size - 1, nullptr),
      loads(2 * n_size, nullptr),
      DGParam_list(2 * n_size)
  {
    buildScaleMicrogridNetwork();
  }

  /**
   * @brief Construct all components of the scaled microgrid network.
   *
   * Builds a microgrid containing @c 2*N_size generators and buses connected
   * in a chain by transmission lines. Loads are connected to every other bus,
   * and a virtual DQ bus is associated with each physical bus.
   *
   * The first two generators use the reference generator parameter set, while
   * all remaining generators use the second parameter set. Line parameters
   * alternate along the network, and the first load uses a different parameter
   * set from the remaining loads.
   *
   * The created components are stored in this network and assigned unique model
   * identifiers using @c model_id_next.
   *
   * @pre @c N_size is greater than zero.
   *
   * @post The network contains @c 2*N_size generators and virtual DQ buses.
   * @post The network contains @c 2*N_size-1 transmission lines connecting
   *       consecutive buses.
   * @post The network contains @c N_size loads connected to every other bus.
   * @post @c DGParam_list contains the parameters associated with each generator.
   * @post @c model_id_next is advanced for every component created.
   *
   * @note Components are dynamically allocated and their pointers are stored
   *       in the corresponding component vectors.
   */
  void buildScaleMicrogridNetwork()
  {
    assert(N_size > 0);
    // Every Bus has the same virtual resistance. This is due to numerical stability as mentioned in the paper.
    ScalarT RN = 1.0e4;

    // TODO: add this as parameters
    // DG Params Vector
    // All DGs have the same set of parameters except for the first two.
    GenParams DG_parms1;
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

    GenParams DG_parms2;
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

    DGParam_list.assign(2 * N_size, DG_parms2);

    // First two generators use parameters 1
    if (DGParam_list.size() >= 1)
    {
      DGParam_list[0] = DG_parms1;
    }
    if (DGParam_list.size() >= 2)
    {
      DGParam_list[1] = DG_parms1;
    }

    // line vector params
    // Every odd line has the same parameters and every even line has the same parameters
    ScalarT              rline1 = 0.23;
    ScalarT              Lline1 = 0.1 / (2.0 * M_PI * 50.0);
    ScalarT              rline2 = 0.35;
    ScalarT              Lline2 = 0.58 / (2.0 * M_PI * 50.0);
    std::vector<ScalarT> rline_list(2 * N_size - 1, 0.0);
    std::vector<ScalarT> Lline_list(2 * N_size - 1, 0.0);
    for (IdxT i = 0; i < rline_list.size(); i++)
    {
      rline_list[i] = (i % 2) ? rline2 : rline1;
      Lline_list[i] = (i % 2) ? Lline2 : Lline1;
    }

    // load parms
    // Only the first load has the same paramaters.
    ScalarT rload1 = 3.0;
    ScalarT Lload1 = 2.0 / (2.0 * M_PI * 50.0);
    ScalarT rload2 = 2.0;
    ScalarT Lload2 = 1.0 / (2.0 * M_PI * 50.0);

    std::vector<ScalarT> rload_list(N_size, rload2);
    std::vector<ScalarT> Lload_list(N_size, Lload2);
    if (rload_list.size() >= 1)
    {
      rload_list[0] = rload1;
      Lload_list[0] = Lload1;
    }

    // Create the reference generator
    auto* dg_ref = new DGGenerator(model_id_next++,
                                   DGParam_list[0],
                                   true,
                                   &dg_signal,
                                   &buses[0]);

    generators[0] = dg_ref;

    // Create the remaining generators.
    for (IdxT i = 1; i < 2 * N_size; i++)
    {
      auto* dg = new DGGenerator(model_id_next++,
                                 DGParam_list[i],
                                 false,
                                 &dg_signal,
                                 &buses[i]);

      generators[i] = dg;
    }

    // Create transmission lines between consecutive buses.
    for (IdxT i = 0; i < 2 * N_size - 1; i++)
    {
      auto* line_model = new Line(model_id_next++,
                                  rline_list[i],
                                  Lline_list[i],
                                  &dg_signal,
                                  &buses[i],
                                  &buses[i + 1]);

      lines[i] = line_model;
    }

    // Create loads on every other bus.
    for (IdxT i = 0; i < N_size; i++)
    {
      auto* load_model = new Load(model_id_next++,
                                  rload_list[i],
                                  Lload_list[i],
                                  &dg_signal,
                                  &buses[2 * i]);

      loads[2 * i] = load_model;
    }

    // Create and Add all the microgrid Virtual DQ Buses
    for (IdxT i = 0; i < 2 * N_size; i++)
    {
      auto* virDQbus_model = new BusDQ(model_id_next++,
                                       RN,
                                       &buses[i]);

      busesDQ[i] = virDQbus_model;
    }
  }
};

#pragma once

#include <array>

#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;

    struct SystemTestData
    {
      using RealT                   = double;
      using IdxT                    = size_t;
      static constexpr RealT SMAX   = 0.6;
      static constexpr RealT VMIN   = 0.95;
      static constexpr RealT VMAX   = 1.05;
      static constexpr RealT LOAD_P = 1.2;
      static constexpr RealT LOAD_Q = 0.3;

      /// Equality rows of the buses in order
      static constexpr std::array<IdxT, 6> BALANCE_ROWS{0, 1, 2, 3, 4, 5};
      /// Inequality rows of the buses in order, then of the rated branch, then
      /// of the unrated branches
      static constexpr std::array<IdxT, 3> MAGNITUDE_ROWS{0, 1, 2};
      static constexpr std::array<IdxT, 2> FLOW_ROWS{3, 4};
      static constexpr std::array<IdxT, 4> UNRATED_ROWS{5, 6, 7, 8};

      /// A point away from the flat start
      static constexpr std::array<RealT, 10> POINT{1.02, 0.0, 0.97, -0.07, 1.01, -0.03, 0.8, 0.2, 0.4, -0.1};

      static OPF::SystemModelData<RealT, IdxT> data()
      {
        OPF::SystemModelData<RealT, IdxT> data;

        for (IdxT number = 1; number <= 3; ++number)
        {
          auto& bus                                     = data.bus.emplace_back();
          bus.number                                    = number;
          bus.parameters[OPF::BusParameters::reference] = number == 1;
          bus.parameters[OPF::BusParameters::Vmin]      = VMIN;
          bus.parameters[OPF::BusParameters::Vmax]      = VMAX;
        }
        data.bus[2].parameters[OPF::BusParameters::B] = 0.1;

        addBranch(data, "branch_1_2", 1, 2, 0.01, 0.1, 0.02);
        addBranch(data, "branch_2_3", 2, 3, 0.01, 0.1, 0.0);
        addBranch(data, "branch_1_3", 1, 3, 0.02, 0.2, 0.0);
        data.branch[0].parameters[OPF::BranchParameters::Smax] = SMAX;

        addGenerator(data, "gen_1", 1, 10.0);
        addGenerator(data, "gen_3", 3, 20.0);

        auto& load                              = data.load.emplace_back();
        load.id                                 = "load_2";
        load.ports[OPF::LoadPorts::bus]         = 2;
        load.parameters[OPF::LoadParameters::P] = LOAD_P;
        load.parameters[OPF::LoadParameters::Q] = LOAD_Q;

        return data;
      }

      static void addBranch(OPF::SystemModelData<RealT, IdxT>& data,
                            const char*                        id,
                            IdxT                               bus1,
                            IdxT                               bus2,
                            RealT                              r,
                            RealT                              x,
                            RealT                              b)
      {
        auto& branch                                = data.branch.emplace_back();
        branch.id                                   = id;
        branch.ports[OPF::BranchPorts::bus1]        = bus1;
        branch.ports[OPF::BranchPorts::bus2]        = bus2;
        branch.parameters[OPF::BranchParameters::R] = r;
        branch.parameters[OPF::BranchParameters::X] = x;
        branch.parameters[OPF::BranchParameters::B] = b;
      }

      static void addGenerator(OPF::SystemModelData<RealT, IdxT>& data, const char* id, IdxT bus, RealT c1)
      {
        auto& generator                                      = data.generator.emplace_back();
        generator.id                                         = id;
        generator.ports[OPF::GeneratorPorts::bus]            = bus;
        generator.parameters[OPF::GeneratorParameters::Pmin] = 0.0;
        generator.parameters[OPF::GeneratorParameters::Pmax] = 2.0;
        generator.parameters[OPF::GeneratorParameters::Qmin] = -1.0;
        generator.parameters[OPF::GeneratorParameters::Qmax] = 1.0;
        generator.parameters[OPF::GeneratorParameters::c1]   = c1;
        generator.parameters[OPF::GeneratorParameters::c2]   = 0.01;
      }
    };
  } // namespace Testing
} // namespace GridKit

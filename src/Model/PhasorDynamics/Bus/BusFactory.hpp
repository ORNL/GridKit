
#pragma once

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarT = double, typename IdxT = int>
    class BusFactory
    {
    public:
      using real_type = typename Model::Evaluator<ScalarT, IdxT>::real_type;
      using BusData   = GridKit::PhasorDynamics::BusData<real_type, IdxT>;
      using BusTypeT  = typename GridKit::PhasorDynamics::BusData<real_type, IdxT>::BusType;

      BusFactory() = delete;

      static BusBase<ScalarT, IdxT>* create(const BusData& data)
      {
        BusBase<ScalarT, IdxT>* bus = nullptr;

        switch (data.bus_type)
        {
        case BusTypeT::DEFAULT:
          bus = new Bus<ScalarT, IdxT>(data);
          break;
        case BusTypeT::SLACK:
          bus = new BusInfinite<ScalarT, IdxT>(data);
          break;
        default:
          // Throw exception
          std::cout << "Bus type " << static_cast<int>(data.bus_type) << " unrecognized.\n";
        }
        return bus;
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit

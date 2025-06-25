
#pragma once

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite/BusInfinite.hpp>
#include <Model/PhasorDynamics/Bus/BusSignal/BusSignal.hpp>

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

      BusFactory() = delete;

      static BusBase<ScalarT, IdxT>* create(const BusData& data)
      {
        BusBase<ScalarT, IdxT>* bus = nullptr;

        switch (data.bus_type)
        {
        case BusData::DEFAULT:
          bus = new Bus<ScalarT, IdxT>(data);
          break;
        case BusData::SLACK:
          bus = new BusInfinite<ScalarT, IdxT>(data);
          break;
        case BusData::SIGNAL:
          bus = new BusSignal<ScalarT, IdxT>(data);
          break;
        default:
          // Throw exception
          std::cout << "Bus type " << data.bus_type << " unrecognized.\n";
        }
        return bus;
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit

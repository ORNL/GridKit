
#pragma once

#include <Model/PhasorDynamics/Bus/BusElectric/BusElectricData.hpp>
#include <Model/PhasorDynamics/Bus/BusElectric/BusNetwork.hpp>
#include <Model/PhasorDynamics/Bus/BusElectric/BusInfinite.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarT = double, typename IdxT = int>
    class BusFactory
    {
    public:
      using real_type = typename BusElectric<ScalarT, IdxT>::real_type;
      using BusElectricData   = GridKit::PhasorDynamics::BusElectricData<real_type, IdxT>;

      BusFactory() = delete;

      static BusElectric<ScalarT, IdxT>* create(const BusElectricData& data)
      {
        BusElectric<ScalarT, IdxT>* bus = nullptr;

        switch (data.bus_type)
        {
        case BusElectricData::DEFAULT:
          bus = new BusNetwork<ScalarT, IdxT>(data);
          break;
        case BusElectricData::SLACK:
          bus = new BusInfinite<ScalarT, IdxT>(data);
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

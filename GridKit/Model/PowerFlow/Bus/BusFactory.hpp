
#pragma once

#include <GridKit/Model/PowerFlow/Bus/BusPQ.hpp>
#include <GridKit/Model/PowerFlow/Bus/BusPV.hpp>
#include <GridKit/Model/PowerFlow/Bus/BusSlack.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{

  template <typename ScalarT = double, typename IdxT = int>
  class BusFactory
  {
  public:
    using RealT   = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BusData = GridKit::PowerFlowData::BusData<RealT, IdxT>;

    BusFactory() = delete;

    static BaseBus<ScalarT, IdxT>* create(BusData& data)
    {
      BaseBus<ScalarT, IdxT>* bus = nullptr;
      switch (data.type)
      {
      case 1:
        bus = new BusPQ<ScalarT, IdxT>(data);
        break;
      case 2:
        bus = new BusPV<ScalarT, IdxT>(data);
        break;
      case 3:
        bus = new BusSlack<ScalarT, IdxT>(data);
        break;
      default:
        // Throw exception
        std::cout << "Bus type " << data.type << " unrecognized.\n";
      }
      return bus;
    }
  };

} // namespace GridKit

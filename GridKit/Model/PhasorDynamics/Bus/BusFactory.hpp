
#pragma once

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageIn.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageOut/BusSignalVoltageOut.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type = double, typename index_type = int>
    class BusFactory
    {
    public:
      using ScalarT        = scalar_type;
      using IdxT           = index_type;
      using RealT          = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using SignalNodeSetT = SignalNodeSet<ScalarT, IdxT>;
      using BusData        = GridKit::PhasorDynamics::BusData<RealT, IdxT>;
      using BusTypeT       = typename GridKit::PhasorDynamics::BusData<RealT, IdxT>::BusType;

      BusFactory() = delete;

      static BusBase<ScalarT, IdxT>* create(const BusData& data, SignalNodeSetT& signal_nodes)
      {
        switch (data.bus_type)
        {
        case BusTypeT::DEFAULT:
          return new Bus<ScalarT, IdxT>(data);
        case BusTypeT::SLACK:
          return new BusInfinite<ScalarT, IdxT>(data);
        case BusTypeT::SIGNAL_VOLTAGE_IN:
        {
          auto* bus = new BusSignalVoltageIn<ScalarT, IdxT>(data);
          bus->getPorts().connect(data, signal_nodes);
          return bus;
        }
        case BusTypeT::SIGNAL_VOLTAGE_OUT:
        {
          auto* bus = new BusSignalVoltageOut<ScalarT, IdxT>(data);
          bus->getPorts().connect(data, signal_nodes);
          return bus;
        }
        default:
          // Throw exception
          ::GridKit::Utilities::Logger::error() << "Bus type " << static_cast<int>(data.bus_type) << " unrecognized.\n";
          return nullptr;
        }
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit


#pragma once

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarT = double, typename IdxT = int>
    class BusFactory
    {
    public:
      using RealT    = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using BusData  = GridKit::PhasorDynamics::BusData<RealT, IdxT>;
      using BusTypeT = typename GridKit::PhasorDynamics::BusData<RealT, IdxT>::BusType;

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
          // TODO: align with BranchFactory::create — throw std::runtime_error on
          // unrecognized bus_type instead of logging and returning nullptr, so
          // misconfiguration fails at construction rather than propagating a
          // null pointer through SystemModel::build.
          std::cout << "Bus type " << static_cast<int>(data.bus_type) << " unrecognized.\n";
        }
        return bus;
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit

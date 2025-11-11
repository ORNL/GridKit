
#pragma once

#include <GridKit/Model/PowerFlow/Bus/BaseBus.hpp>
#include <GridKit/Model/PowerFlow/Generator/GeneratorPQ.hpp>
#include <GridKit/Model/PowerFlow/Generator/GeneratorPV.hpp>
#include <GridKit/Model/PowerFlow/Generator/GeneratorSlack.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{

  template <typename ScalarT = double, typename IdxT = int>
  class GeneratorFactory
  {
  public:
    using RealT   = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using GenData = GridKit::PowerFlowData::GenData<RealT, IdxT>;

    GeneratorFactory() = delete;

    static GeneratorBase<ScalarT, IdxT>* create(BaseBus<ScalarT, IdxT>* bus, GenData& data)
    {
      GeneratorBase<ScalarT, IdxT>* gen = nullptr;
      switch (bus->BusType())
      {
      case 1:
        gen = new GeneratorPQ<ScalarT, IdxT>(bus, data);
        break;
      case 2:
        gen = new GeneratorPV<ScalarT, IdxT>(bus, data);
        break;
      case 3:
        gen = new GeneratorSlack<ScalarT, IdxT>(bus, data);
        break;
      default:
        // Throw exception
        std::cout << "Generator type " << bus->BusType() << " unrecognized.\n";
      }
      return gen;
    }
  };

} // namespace GridKit

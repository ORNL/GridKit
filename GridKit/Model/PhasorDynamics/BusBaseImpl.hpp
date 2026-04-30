#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElementImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename ScalarP, typename IdxP>
    BusBase<ScalarP, IdxP>::BusBase(const BusData<RealT, IdxT>& data)
      : bus_id_(data.bus_id)
    {
    }

    template <typename ScalarP, typename IdxP>
    BusBase<ScalarP, IdxP>::~BusBase()
    {
    }
  } // namespace PhasorDynamics
} // namespace GridKit

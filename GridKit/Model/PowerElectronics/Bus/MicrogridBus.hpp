#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class MicrogridBus : public NodeBase<ScalarT, IdxT>
    {
    public:
      MicrogridBus() : NodeBase<ScalarT, IdxT>(2, 0)
      {
      }
    };
  } // namespace PowerElectronics
} // namespace GridKit
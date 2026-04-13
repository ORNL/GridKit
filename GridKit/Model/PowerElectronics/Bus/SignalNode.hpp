#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class SignalNode : public NodeBase<ScalarT, IdxT>
    {
    public:
      SignalNode()
        : NodeBase<ScalarT, IdxT>(1, 0)
      {
      }
    };
  } // namespace PowerElectronics
} // namespace GridKit

#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class DGSignal : public NodeBase<ScalarT, IdxT>
    {
    public:
      DGSignal() : NodeBase<ScalarT, IdxT>(1, 0)
      {
      }
    };
  } // namespace PowerElectronics
} // namespace GridKit
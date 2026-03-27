#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class Bus : public NodeBase<ScalarT, IdxT>
    {
    public:
      Bus() : NodeBase<ScalarT, IdxT>(1, 0)
      {
      }
    };
  } // namespace PowerElectronics
} // namespace GridKit
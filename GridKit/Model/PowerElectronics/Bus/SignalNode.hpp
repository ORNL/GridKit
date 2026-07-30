#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename scalar_type, typename index_type>
    class SignalNode : public NodeBase<scalar_type, index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;

      SignalNode()
        : NodeBase<ScalarT, IdxT>(1, 0)
      {
      }
    };
  } // namespace PowerElectronics
} // namespace GridKit

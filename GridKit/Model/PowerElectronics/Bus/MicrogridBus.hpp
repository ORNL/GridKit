#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename scalar_type, typename index_type>
    class MicrogridBus : public NodeBase<scalar_type, index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;

      MicrogridBus()
        : NodeBase<ScalarT, IdxT>(2, 0)
      {
      }
    };
  } // namespace PowerElectronics
} // namespace GridKit

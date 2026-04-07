#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class InfiniteBus : public NodeBase<ScalarT, IdxT>
    {
      using NodeBase<ScalarT, IdxT>::y;

    public:
      InfiniteBus(ScalarT voltage) : NodeBase<ScalarT, IdxT>(0, 1), voltage_(voltage)
      {
      }

      int initialize() final
      {
        if (int err_code = NodeBase<ScalarT, IdxT>::initialize())
          return err_code;

        y()[0] = voltage_;

        return 0;
      }

    private:
      ScalarT voltage_;
    };
  } // namespace PowerElectronics
} // namespace GridKit
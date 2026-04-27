#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class GroundedBus : public NodeBase<ScalarT, IdxT>
    {
      using ExternalConnection = typename CircuitComponent<ScalarT, IdxT>::ExternalConnection;

      using NodeBase<ScalarT, IdxT>::y_ext_;
      using NodeBase<ScalarT, IdxT>::yp_ext_;
      using NodeBase<ScalarT, IdxT>::f_ext_;

    public:
      GroundedBus(ScalarT voltage)
        : NodeBase<ScalarT, IdxT>(0, 1), voltage_(voltage)
      {
      }

      int allocate() final
      {
        if (int err_code = NodeBase<ScalarT, IdxT>::allocate())
          return err_code;

        this->setExternalConnectionNodes(0, ExternalConnection{.y_ = &voltage_, .yp_ = &dummy_, .f_ = &dummy_, .idx_ = INVALID_INDEX<IdxT>});

        return 0;
      }

    private:
      ScalarT voltage_;
      ScalarT dummy_;
    };
  } // namespace PowerElectronics
} // namespace GridKit

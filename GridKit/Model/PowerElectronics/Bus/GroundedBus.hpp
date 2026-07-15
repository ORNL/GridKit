#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class GroundedBus : public NodeBase<ScalarT, IdxT>
    {
      using NodeBase<ScalarT, IdxT>::y_ext_;
      using NodeBase<ScalarT, IdxT>::yp_ext_;
      using NodeBase<ScalarT, IdxT>::f_ext_;

    public:
      GroundedBus(ScalarT voltage)
        : NodeBase<ScalarT, IdxT>(0, 1), voltage_(voltage)
      {
      }

      int initialize() final
      {
        if (int err_code = NodeBase<ScalarT, IdxT>::initialize())
          return err_code;

        auto* y_data = y().getData();
        y_data[0]    = voltage_;
        y().setDataUpdated();

        return 0;
      }

      int allocate() final
      {
        if (int err_code = NodeBase<ScalarT, IdxT>::allocate())
          return err_code;

        this->setExternalConnectionNodes(0, ExternalConnection<ScalarT, IdxT>{.y_ = &voltage_, .yp_ = &dummy_, .f_ = &dummy_, .idx_ = INVALID_INDEX<IdxT>});

        return 0;
      }

    private:
      ScalarT voltage_;
      ScalarT dummy_ = 0;
    };
  } // namespace PowerElectronics
} // namespace GridKit

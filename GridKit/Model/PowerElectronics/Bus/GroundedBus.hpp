#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class GroundedBus : public NodeBase<ScalarT, IdxT>
    {
      using NodeBase<ScalarT, IdxT>::y;

    public:
      GroundedBus(ScalarT voltage)
        : NodeBase<ScalarT, IdxT>(0, 1), voltage_(voltage)
      {
      }

      int initialize() final
      {
        if (int err_code = NodeBase<ScalarT, IdxT>::initialize())
          return err_code;

        y()[0] = voltage_;

        return 0;
      }

      int allocate() final
      {
        if (int err_code = NodeBase<ScalarT, IdxT>::allocate())
          return err_code;

        this->setExternalConnectionNodes(0, INVALID_INDEX<IdxT>);

        return 0;
      }

    private:
      ScalarT voltage_;
    };
  } // namespace PowerElectronics
} // namespace GridKit

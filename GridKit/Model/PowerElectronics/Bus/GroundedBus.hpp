#pragma once

#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename scalar_type, typename index_type>
    class GroundedBus : public NodeBase<scalar_type, index_type>
    {
      using NodeBase<scalar_type, index_type>::y;

    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;

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

        this->setExternalConnectionNodes(0, INVALID_INDEX<IdxT>);

        return 0;
      }

    private:
      ScalarT voltage_;
    };
  } // namespace PowerElectronics
} // namespace GridKit

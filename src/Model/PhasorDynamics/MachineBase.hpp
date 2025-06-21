#pragma once

#include <vector>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief MachineBase model implementation base class.
     *
     */
    template <class ScalarT, typename IdxT>
    class MachineBase
    {
    public:
      virtual ScalarT speed() = 0;
      virtual ScalarT get_torque() = 0;
    };

  } // namespace PhasorDynamics
} // namespace GridKit

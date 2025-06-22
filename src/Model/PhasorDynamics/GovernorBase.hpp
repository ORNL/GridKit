#pragma once

#include <vector>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief GovernorBase model implementation base class.
     *
     */
    template <class ScalarT, typename IdxT>
    class GovernorBase
    {
    public:
      virtual ScalarT& Pmech() = 0;
    };

  } // namespace PhasorDynamics
} // namespace GridKit

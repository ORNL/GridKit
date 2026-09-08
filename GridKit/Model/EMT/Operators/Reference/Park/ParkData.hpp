#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class ParkParameters
    {
      inverse, ///< Apply the inverse Park transformation
    };

    enum class ParkInputs : size_t
    {
      u1,    ///< \f$u_1\f$ Input component
      u2,    ///< \f$u_2\f$ Input component
      u3,    ///< \f$u_3\f$ Input component
      theta, ///< \f$\theta\f$ Electrical angle [rad]
      SIZE,
    };

    enum class ParkOutputs : size_t
    {
      y1, ///< \f$y_1\f$ Output component
      y2, ///< \f$y_2\f$ Output component
      y3, ///< \f$y_3\f$ Output component
      SIZE,
    };

    enum class ParkMonitorableVariables
    {
      y1, ///< \f$y_1\f$ Output component
      y2, ///< \f$y_2\f$ Output component
      y3, ///< \f$y_3\f$ Output component
    };

    template <typename real_type, typename index_type>
    struct ParkData : public ComponentData<real_type,
                                           index_type,
                                           ParkParameters,
                                           ParkInputs,
                                           ParkOutputs,
                                           ParkMonitorableVariables>
    {
      ParkData() = default;

      using Parameters           = ParkParameters;
      using Inputs               = ParkInputs;
      using Outputs              = ParkOutputs;
      using MonitorableVariables = ParkMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit

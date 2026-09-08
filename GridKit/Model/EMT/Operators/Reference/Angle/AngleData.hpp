#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class AngleParameters
    {
    };

    enum class AngleInputs : size_t
    {
      omega, ///< \f$\omega\f$ Electrical angular frequency [rad/s]
      SIZE,
    };

    enum class AngleOutputs : size_t
    {
      theta, ///< \f$\theta\f$ Electrical angle [rad]
      SIZE,
    };

    enum class AngleMonitorableVariables
    {
      theta, ///< \f$\theta\f$ Electrical angle [rad]
    };

    template <typename real_type, typename index_type>
    struct AngleData : public ComponentData<real_type,
                                            index_type,
                                            AngleParameters,
                                            AngleInputs,
                                            AngleOutputs,
                                            AngleMonitorableVariables>
    {
      AngleData() = default;

      using Parameters           = AngleParameters;
      using Inputs               = AngleInputs;
      using Outputs              = AngleOutputs;
      using MonitorableVariables = AngleMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit

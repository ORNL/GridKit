#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class OuterVoltageControlParameters
      {
        C,
        Kp,
        Ki,
        Kaw,
      };

      enum class OuterVoltageControlInputs : size_t
      {
        vrefd,
        vrefq,
        vd,
        vq,
        igd,
        igq,
        omega,
        ilimd,
        ilimq,
        SIZE,
      };

      enum class OuterVoltageControlOutputs : size_t
      {
        irefd,
        irefq,
        SIZE,
      };

      enum class OuterVoltageControlMonitorableVariables
      {
        etad,
        etaq,
        irefd,
        irefq,
      };

      template <typename real_type, typename index_type>
      struct OuterVoltageControlData : public ComponentData<real_type,
                                                            index_type,
                                                            OuterVoltageControlParameters,
                                                            OuterVoltageControlInputs,
                                                            OuterVoltageControlOutputs,
                                                            OuterVoltageControlMonitorableVariables>
      {
        OuterVoltageControlData() = default;

        using Parameters           = OuterVoltageControlParameters;
        using Inputs               = OuterVoltageControlInputs;
        using Outputs              = OuterVoltageControlOutputs;
        using MonitorableVariables = OuterVoltageControlMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

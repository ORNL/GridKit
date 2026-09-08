#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class InnerCurrentControlParameters
      {
        L,
        Kp,
        Ki,
        Kaw,
        Imax,
        Mmax,
      };

      enum class InnerCurrentControlInputs : size_t
      {
        vd,
        vq,
        id,
        iq,
        irefd,
        irefq,
        omega,
        vdc,
        SIZE,
      };

      enum class InnerCurrentControlOutputs : size_t
      {
        ilimd,
        ilimq,
        ud,
        uq,
        SIZE,
      };

      enum class InnerCurrentControlMonitorableVariables
      {
        xid,
        xiq,
        ilimd,
        ilimq,
        ud,
        uq,
      };

      template <typename real_type, typename index_type>
      struct InnerCurrentControlData : public ComponentData<real_type,
                                                            index_type,
                                                            InnerCurrentControlParameters,
                                                            InnerCurrentControlInputs,
                                                            InnerCurrentControlOutputs,
                                                            InnerCurrentControlMonitorableVariables>
      {
        InnerCurrentControlData() = default;

        using Parameters           = InnerCurrentControlParameters;
        using Inputs               = InnerCurrentControlInputs;
        using Outputs              = InnerCurrentControlOutputs;
        using MonitorableVariables = InnerCurrentControlMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

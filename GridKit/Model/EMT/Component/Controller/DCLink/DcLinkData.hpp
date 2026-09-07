#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class DcLinkParameters
      {
        C,
      };

      enum class DcLinkInputs : size_t
      {
        isrc,
        idc,
        SIZE,
      };

      enum class DcLinkOutputs : size_t
      {
        vdc,
        SIZE,
      };

      enum class DcLinkMonitorableVariables
      {
        vdc,
        isrc,
        idc,
        energy,
      };

      template <typename real_type, typename index_type>
      struct DcLinkData : public ComponentData<real_type,
                                               index_type,
                                               DcLinkParameters,
                                               DcLinkInputs,
                                               DcLinkOutputs,
                                               DcLinkMonitorableVariables>
      {
        DcLinkData() = default;

        using Parameters           = DcLinkParameters;
        using Inputs               = DcLinkInputs;
        using Outputs              = DcLinkOutputs;
        using MonitorableVariables = DcLinkMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

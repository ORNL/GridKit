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
        C, ///< \f$C\f$ Capacitance [F]
      };

      enum class DcLinkInputs : size_t
      {
        isrc, ///< \f$i_{\mathrm{src}}\f$ Source current [A]
        idc,  ///< \f$i_{\mathrm{dc}}\f$ Converter current [A]
        SIZE,
      };

      enum class DcLinkOutputs : size_t
      {
        vdc, ///< \f$v_{\mathrm{dc}}\f$ DC-link voltage [V]
        SIZE,
      };

      enum class DcLinkMonitorableVariables
      {
        vdc,    ///< \f$v_{\mathrm{dc}}\f$ DC-link voltage [V]
        isrc,   ///< \f$i_{\mathrm{src}}\f$ Source current [A]
        idc,    ///< \f$i_{\mathrm{dc}}\f$ Converter current [A]
        energy, ///< \f$E\f$ Capacitor energy [J]
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

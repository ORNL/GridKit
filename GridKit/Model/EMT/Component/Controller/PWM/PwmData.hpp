#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class PwmParameters
      {
        M,         ///< \f$M\f$ Prescribed modulation index [-]
        fm,        ///< \f$f_{\mathrm{m}}\f$ Modulation frequency [Hz]
        fc,        ///< \f$f_{\mathrm{c}}\f$ Carrier frequency [Hz]
        alignment, ///< \f$\alpha\f$ Pulse alignment [-]
      };

      enum class PwmInputs : size_t
      {
        ma, ///< \f$m_a\f$ Phase modulation command [-]
        mb, ///< \f$m_b\f$ Phase modulation command [-]
        mc, ///< \f$m_c\f$ Phase modulation command [-]
        SIZE,
      };

      enum class PwmOutputs : size_t
      {
        sa, ///< \f$s_a\f$ Phase switching function [-]
        sb, ///< \f$s_b\f$ Phase switching function [-]
        sc, ///< \f$s_c\f$ Phase switching function [-]
        SIZE,
      };

      enum class PwmMonitorableVariables
      {
        sa, ///< \f$s_a\f$ Phase switching function [-]
        sb, ///< \f$s_b\f$ Phase switching function [-]
        sc, ///< \f$s_c\f$ Phase switching function [-]
      };

      template <typename real_type, typename index_type>
      struct PwmData : public ComponentData<real_type,
                                            index_type,
                                            PwmParameters,
                                            PwmInputs,
                                            PwmOutputs,
                                            PwmMonitorableVariables>
      {
        PwmData() = default;

        using Parameters           = PwmParameters;
        using Inputs               = PwmInputs;
        using Outputs              = PwmOutputs;
        using MonitorableVariables = PwmMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

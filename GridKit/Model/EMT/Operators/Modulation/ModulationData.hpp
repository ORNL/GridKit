#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class ModulationParameters
    {
      Mmax, ///< \f$M^{\max}\f$ Sinusoidal modulation limit [-]
    };

    enum class ModulationInputs : size_t
    {
      ud,  ///< \f$u_d\f$ Converter voltage command [V]
      uq,  ///< \f$u_q\f$ Converter voltage command [V]
      vdc, ///< \f$v_{\mathrm{dc}}\f$ DC-link voltage [V]
      SIZE,
    };

    enum class ModulationOutputs : size_t
    {
      md,    ///< \f$m_d\f$ Modulation command [-]
      mq,    ///< \f$m_q\f$ Modulation command [-]
      ulimd, ///< \f$u_d^{\mathrm{lim}}\f$ Limited voltage command [V]
      ulimq, ///< \f$u_q^{\mathrm{lim}}\f$ Limited voltage command [V]
      SIZE,
    };

    enum class ModulationMonitorableVariables
    {
      md,    ///< \f$m_d\f$ Modulation command [-]
      mq,    ///< \f$m_q\f$ Modulation command [-]
      ulimd, ///< \f$u_d^{\mathrm{lim}}\f$ Limited voltage command [V]
      ulimq, ///< \f$u_q^{\mathrm{lim}}\f$ Limited voltage command [V]
    };

    template <typename real_type, typename index_type>
    struct ModulationData : public ComponentData<real_type,
                                                 index_type,
                                                 ModulationParameters,
                                                 ModulationInputs,
                                                 ModulationOutputs,
                                                 ModulationMonitorableVariables>
    {
      ModulationData() = default;

      using Parameters           = ModulationParameters;
      using Inputs               = ModulationInputs;
      using Outputs              = ModulationOutputs;
      using MonitorableVariables = ModulationMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit

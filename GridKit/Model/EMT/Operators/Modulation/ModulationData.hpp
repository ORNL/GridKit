#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class ModulationParameters
    {
    };

    enum class ModulationInputs : size_t
    {
      ua,  ///< \f$u_a\f$ Phase voltage command [V]
      ub,  ///< \f$u_b\f$ Phase voltage command [V]
      uc,  ///< \f$u_c\f$ Phase voltage command [V]
      vdc, ///< \f$v_{\mathrm{dc}}\f$ DC-link voltage [V]
      SIZE,
    };

    enum class ModulationOutputs : size_t
    {
      ma, ///< \f$m_a\f$ Phase modulation command [-]
      mb, ///< \f$m_b\f$ Phase modulation command [-]
      mc, ///< \f$m_c\f$ Phase modulation command [-]
      SIZE,
    };

    enum class ModulationMonitorableVariables
    {
      ma, ///< \f$m_a\f$ Phase modulation command [-]
      mb, ///< \f$m_b\f$ Phase modulation command [-]
      mc, ///< \f$m_c\f$ Phase modulation command [-]
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

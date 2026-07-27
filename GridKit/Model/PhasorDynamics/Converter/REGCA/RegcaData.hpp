/**
 * @file RegcaData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the REGCA converter model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Parameter keys for the REGCA converter model.
      enum class RegcaParameters
      {
        p0,     ///< Initial active power injection on system base
        q0,     ///< Initial reactive power injection on system base
        mva,    ///< MVA base of the REGCA model
        Tg,     ///< Converter current-control lag time constant
        TM,     ///< Terminal voltage sensor time constant
        Rqmax,  ///< Reactive-current recovery positive rate limit
        Rqmin,  ///< Reactive-current recovery negative rate limit
        Rpmax,  ///< Active-current magnitude recovery rate limit
        sL,     ///< LVPL switch
        IL1,    ///< LVPL upper-current ceiling
        VL0,    ///< LVPL zero-crossing voltage
        VL1,    ///< LVPL upper breakpoint voltage
        VA0,    ///< LVACM lower breakpoint voltage
        VA1,    ///< LVACM upper breakpoint voltage
        Vhvmax, ///< Terminal-voltage ceiling for HV reactive management

        // Optional PowerWorld compatibility fields, accepted and unused.
        Qmin, ///< Unused compatibility field
        Khv,  ///< Unused compatibility field
        Xe    ///< Unused compatibility field
      };

      /// Buses for the REGCA converter model.
      enum class RegcaBuses : size_t
      {
        bus, ///< Terminal bus ID
        SIZE
      };

      /// Signal inputs for the REGCA converter model.
      enum class RegcaSignalInputs : size_t
      {
        ipcmd, ///< Optional active-current command signal ID
        iqcmd, ///< Optional reactive-current command signal ID
        SIZE
      };

      /// Signal outputs for the REGCA converter model.
      enum class RegcaSignalOutputs : size_t
      {
        ibranchr, ///< Optional branch-current real-component output signal ID
        ibranchi, ///< Optional branch-current imaginary-component output signal ID
        pbranch,  ///< Optional branch active-power output signal ID
        qbranch,  ///< Optional branch reactive-power output signal ID
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class RegcaMonitorableVariables
      {
        ir, ///< Branch-current real component on system base
        ii, ///< Branch-current imaginary component on system base
        p,  ///< Branch active power on system base
        q   ///< Branch reactive power on system base
      };

      template <typename real_type, typename index_type>
      struct RegcaData : public ComponentData<real_type,
                                              index_type,
                                              RegcaParameters,
                                              RegcaBuses,
                                              RegcaSignalInputs,
                                              RegcaSignalOutputs,
                                              RegcaMonitorableVariables>
      {
        RegcaData() = default;

        using Parameters           = RegcaParameters;
        using Buses                = RegcaBuses;
        using SignalInputs         = RegcaSignalInputs;
        using SignalOutputs        = RegcaSignalOutputs;
        using MonitorableVariables = RegcaMonitorableVariables;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

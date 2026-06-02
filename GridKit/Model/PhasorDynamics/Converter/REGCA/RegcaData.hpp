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
        P0,    ///< Initial active power on system base
        Q0,    ///< Initial reactive power on system base
        mva,   ///< MVA base of the REGCA model
        Tg,    ///< Converter current-control lag time constant
        TM,    ///< Terminal voltage sensor time constant
        Rqmax, ///< Reactive-current recovery positive rate limit
        Rqmin, ///< Reactive-current recovery negative rate limit
        Rpmax, ///< Active-current magnitude recovery rate limit
        sL,    ///< LVPL switch
        IL1,   ///< LVPL upper-current ceiling
        VL0,   ///< LVPL zero-crossing voltage
        VL1,   ///< LVPL upper breakpoint voltage
        VA0,   ///< LVACM lower breakpoint voltage
        VA1,   ///< LVACM upper breakpoint voltage
        Vhvmax ///< Terminal-voltage ceiling for HV reactive management
      };

      /// Ports for the REGCA converter model.
      enum class RegcaPorts
      {
        bus,      ///< Terminal bus ID
        ipcmd,    ///< Optional active-current command signal ID
        iqcmd,    ///< Optional reactive-current command signal ID
        ibranchr, ///< Optional real current measurement output signal ID
        ibranchi, ///< Optional imaginary current measurement output signal ID
        pbranch,  ///< Optional active-power measurement output signal ID
        qbranch   ///< Optional reactive-power measurement output signal ID
      };

      /// Variables available through the monitor interface.
      enum class RegcaMonitorableVariables
      {
        ir,      ///< Real current injection on converter base
        ii,      ///< Imaginary current injection on converter base
        p,       ///< Active power injection on system base
        q,       ///< Reactive power injection on system base
        vt,      ///< Terminal voltage magnitude
        vm,      ///< Filtered terminal voltage
        ip,      ///< Active-current state
        iq,      ///< Reactive-current state
        iqextra, ///< HVRCM extra reactive current
        il,      ///< LVPL upper-limit current curve
        lp,      ///< Active-current lower rate bound
        up       ///< Active-current upper rate bound
      };

      template <typename real_type, typename index_type>
      struct RegcaData : public ComponentData<real_type,
                                              index_type,
                                              RegcaParameters,
                                              RegcaPorts,
                                              RegcaMonitorableVariables>
      {
        RegcaData() = default;

        using Parameters           = RegcaParameters;
        using Ports                = RegcaPorts;
        using MonitorableVariables = RegcaMonitorableVariables;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

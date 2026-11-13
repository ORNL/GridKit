/**
 * @file RegcaData.hpp
 * @brief Modeling data for optimal power flow REGCA converters.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Converter
    {
      /// Parameters for a REGCA converter, those its steady state reads
      enum class RegcaParameters : size_t
      {
        mva,    ///< \f$S^\mathrm{base}\f$ REGCA component power base [MVA]
        sL,     ///< \f$s_L\f$ LVPL switch
        IL1,    ///< \f$I_{L1}\f$ LVPL upper-current ceiling [p.u.]
        VL1,    ///< \f$V_{L1}\f$ LVPL upper breakpoint voltage [p.u.]
        VA1,    ///< \f$V_{A1}\f$ LVACM upper breakpoint voltage [p.u.]
        Vhvmax, ///< \f$V_\mathrm{hv}^{\max}\f$ Terminal-voltage ceiling [p.u.]
      };

      /// Ports of a REGCA converter
      enum class RegcaPorts : size_t
      {
        bus, ///< Port bus
      };

      /// Signals of a REGCA converter
      enum class RegcaSignals : size_t
      {
        ipcmd,   ///< Active-current command
        iqcmd,   ///< Reactive-current command
        pbranch, ///< Branch active power
        qbranch, ///< Branch reactive power
      };

      /// Data of a `Regca`
      template <typename real_type, typename index_type>
      using RegcaData = ComponentData<real_type, index_type, RegcaParameters, RegcaPorts, RegcaSignals>;
    } // namespace Converter
  } // namespace OptimalPowerFlow
} // namespace GridKit

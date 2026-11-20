/**
 * @file LVRTData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT low-voltage ride-through.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of an `LVRT`, each of size 1 and in per unit
    enum class LVRTInputs : size_t
    {
      vt,    ///< \f$v_t\f$ Terminal voltage
      ipcmd, ///< \f$i_p^\mathrm{cmd}\f$ Active-current order
    };

    /// Parameters of an `LVRT`
    enum class LVRTParameters : size_t
    {
      TM,    ///< \f$T_M\f$ Time constant of the current ceiling [s]
      VL0,   ///< \f$V_{L0}\f$ Voltage at which the ceiling target is zero [p.u.]
      VL1,   ///< \f$V_{L1}\f$ Voltage at which the ceiling target is full [p.u.]
      IL1,   ///< \f$I_{L1}\f$ Full ceiling target [p.u.]
      Rpmax, ///< \f$R_p^{\max}\f$ Maximum recovery rate of the ceiling [p.u./s]
      VA0,   ///< \f$V_{A0}\f$ Voltage at which no current is available [p.u.]
      VA1,   ///< \f$V_{A1}\f$ Voltage at which the whole ceiling is available [p.u.]
    };

    /// Data of an `LVRT`
    template <typename RealT, typename IdxT>
    struct LVRTData
    {
      using Parameters = LVRTParameters;

      RealT TM{ZERO<RealT>};
      RealT VL0{ZERO<RealT>};
      RealT VL1{ZERO<RealT>};
      RealT IL1{ZERO<RealT>};
      RealT Rpmax{ZERO<RealT>};
      RealT VA0{ZERO<RealT>};
      RealT VA1{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::TM, &LVRTData::TM},
                                             Parameter{Parameters::VL0, &LVRTData::VL0},
                                             Parameter{Parameters::VL1, &LVRTData::VL1},
                                             Parameter{Parameters::IL1, &LVRTData::IL1},
                                             Parameter{Parameters::Rpmax, &LVRTData::Rpmax},
                                             Parameter{Parameters::VA0, &LVRTData::VA0},
                                             Parameter{Parameters::VA1, &LVRTData::VA1}};
    };
  } // namespace EMT
} // namespace GridKit

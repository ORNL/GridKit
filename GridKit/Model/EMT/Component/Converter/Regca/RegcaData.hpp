/**
 * @file RegcaData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT renewable energy generator/converter models.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Terminals of a `Regca`, each of size 3
    enum class RegcaTerminals : size_t
    {
      bus, ///< Bus of the converter
    };

    /// Inputs of a `Regca`, each of size 1
    enum class RegcaInputs : size_t
    {
      vt,    ///< \f$V_t\f$ Terminal voltage, line-to-line RMS [V]
      theta, ///< \f$\theta\f$ Angle of the terminal voltage [rad]
      ipcmd, ///< \f$i_p^\mathrm{cmd}\f$ Active-current order [p.u.]
      iqcmd, ///< \f$i_q^\mathrm{cmd}\f$ Reactive-current order [p.u.]
    };

    /// Parameters of a `Regca`
    enum class RegcaParameters : size_t
    {
      V,      ///< \f$V\f$ Rated line-to-line RMS voltage [V]
      mva,    ///< \f$S\f$ Three-phase device rating [MVA]
      Tg,     ///< \f$T_g\f$ Time constant of the currents [s]
      TM,     ///< \f$T_M\f$ Time constant of the current ceiling [s]
      Rqmax,  ///< \f$R_q^{\max}\f$ Upper rate limit of the reactive current, none when zero [p.u./s]
      Rqmin,  ///< \f$R_q^{\min}\f$ Lower rate limit of the reactive current, none when zero [p.u./s]
      Rpmax,  ///< \f$R_p^{\max}\f$ Maximum recovery rate of the ceiling [p.u./s]
      IL1,    ///< \f$I_{L1}\f$ Full ceiling target [p.u.]
      VL0,    ///< \f$V_{L0}\f$ Voltage at which the ceiling target is zero [p.u.]
      VL1,    ///< \f$V_{L1}\f$ Voltage at which the ceiling target is full [p.u.]
      VA0,    ///< \f$V_{A0}\f$ Voltage at which no current is available [p.u.]
      VA1,    ///< \f$V_{A1}\f$ Voltage at which the whole ceiling is available [p.u.]
      Vhvmax, ///< \f$V_\mathrm{hv}^{\max}\f$ Voltage above which reactive current is absorbed [p.u.]
      Khv,    ///< \f$K_\mathrm{hv}\f$ Gain of the voltage response [p.u.]
      Qmin,   ///< \f$Q^{\min}\f$ Lower bound of the reduced order [p.u.]
    };

    /// Data of a `Regca`
    template <typename RealT, typename IdxT>
    struct RegcaData
    {
      using Parameters = RegcaParameters;

      RealT V{ZERO<RealT>};
      RealT mva{ZERO<RealT>};
      RealT Tg{ZERO<RealT>};
      RealT TM{ZERO<RealT>};
      RealT Rqmax{ZERO<RealT>};
      RealT Rqmin{ZERO<RealT>};
      RealT Rpmax{ZERO<RealT>};
      RealT IL1{ZERO<RealT>};
      RealT VL0{ZERO<RealT>};
      RealT VL1{ZERO<RealT>};
      RealT VA0{ZERO<RealT>};
      RealT VA1{ZERO<RealT>};
      RealT Vhvmax{ZERO<RealT>};
      RealT Khv{ZERO<RealT>};
      RealT Qmin{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::V, &RegcaData::V},
                                             Parameter{Parameters::mva, &RegcaData::mva},
                                             Parameter{Parameters::Tg, &RegcaData::Tg},
                                             Parameter{Parameters::TM, &RegcaData::TM},
                                             Parameter{Parameters::Rqmax, &RegcaData::Rqmax},
                                             Parameter{Parameters::Rqmin, &RegcaData::Rqmin},
                                             Parameter{Parameters::Rpmax, &RegcaData::Rpmax},
                                             Parameter{Parameters::IL1, &RegcaData::IL1},
                                             Parameter{Parameters::VL0, &RegcaData::VL0},
                                             Parameter{Parameters::VL1, &RegcaData::VL1},
                                             Parameter{Parameters::VA0, &RegcaData::VA0},
                                             Parameter{Parameters::VA1, &RegcaData::VA1},
                                             Parameter{Parameters::Vhvmax, &RegcaData::Vhvmax},
                                             Parameter{Parameters::Khv, &RegcaData::Khv},
                                             Parameter{Parameters::Qmin, &RegcaData::Qmin}};
    };
  } // namespace EMT
} // namespace GridKit

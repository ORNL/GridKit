/**
 * @file HVRCMData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT high-voltage reactive current management.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of an `HVRCM`, each of size 1 and in per unit
    enum class HVRCMInputs : size_t
    {
      vt,    ///< \f$v_t\f$ Terminal voltage
      iqcmd, ///< \f$i_q^\mathrm{cmd}\f$ Reactive-current order
    };

    /// Parameters of an `HVRCM`
    enum class HVRCMParameters : size_t
    {
      Vhvmax, ///< \f$V_\mathrm{hv}^{\max}\f$ Voltage above which reactive current is absorbed [p.u.]
      Khv,    ///< \f$K_\mathrm{hv}\f$ Gain of the voltage response [p.u.]
      Qmin,   ///< \f$Q^{\min}\f$ Lower bound of the reduced order [p.u.]
    };

    /// Data of an `HVRCM`
    template <typename RealT, typename IdxT>
    struct HVRCMData
    {
      using Parameters = HVRCMParameters;

      RealT Vhvmax{ZERO<RealT>};
      RealT Khv{ZERO<RealT>};
      RealT Qmin{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::Vhvmax, &HVRCMData::Vhvmax},
                                             Parameter{Parameters::Khv, &HVRCMData::Khv},
                                             Parameter{Parameters::Qmin, &HVRCMData::Qmin}};
    };
  } // namespace EMT
} // namespace GridKit

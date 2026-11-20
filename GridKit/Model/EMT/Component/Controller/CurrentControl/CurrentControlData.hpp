/**
 * @file CurrentControlData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT current controls.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `CurrentControl`
    enum class CurrentControlInputs : size_t
    {
      vpos,  ///< \f$v_+\f$ Positive DC terminal voltage [V]
      vneg,  ///< \f$v_-\f$ Negative DC terminal voltage [V]
      vdq,   ///< \f$\mathbf{v}_\mathrm{dq}\f$ Voltages, d then q [V]
      idq,   ///< \f$\mathbf{i}_\mathrm{dq}\f$ Currents, d then q [A]
      omega, ///< \f$\omega\f$ Frequency [rad/s]
      axis,  ///< \f$(c, s)\f$ Active-current unit direction in dq coordinates
      ip,    ///< \f$i_p\f$ Active-current order [p.u.]
      iq,    ///< \f$i_q\f$ Reactive-current order [p.u.]
    };

    /// Parameters of a `CurrentControl`
    enum class CurrentControlParameters : size_t
    {
      V,    ///< \f$V\f$ Rated line-to-line RMS voltage [V]
      mva,  ///< \f$S\f$ Three-phase device rating [MVA]
      Kp,   ///< \f$K_P\f$ Proportional gain [Ohm]
      Ki,   ///< \f$K_I\f$ Integral gain [Ohm/s]
      L,    ///< \f$L\f$ Decoupling inductance [H]
      Mmax, ///< \f$M^{\max}\f$ Maximum modulation magnitude
      C,    ///< \f$C\f$ Filter capacitance for grid-current commands [F]
    };

    /// Data of a `CurrentControl`
    template <typename RealT, typename IdxT>
    struct CurrentControlData
    {
      using Parameters = CurrentControlParameters;

      RealT V{ZERO<RealT>};
      RealT mva{ZERO<RealT>};
      RealT Kp{ZERO<RealT>};
      RealT Ki{ZERO<RealT>};
      RealT L{ZERO<RealT>};
      RealT Mmax{ONE<RealT>};
      RealT C{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::V, &CurrentControlData::V},
                                             Parameter{Parameters::mva, &CurrentControlData::mva},
                                             Parameter{Parameters::Kp, &CurrentControlData::Kp},
                                             Parameter{Parameters::Ki, &CurrentControlData::Ki},
                                             Parameter{Parameters::L, &CurrentControlData::L},
                                             Parameter{Parameters::Mmax, &CurrentControlData::Mmax},
                                             Parameter{Parameters::C, &CurrentControlData::C}};
    };
  } // namespace EMT
} // namespace GridKit

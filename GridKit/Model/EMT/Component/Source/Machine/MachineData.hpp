/**
 * @file MachineData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT synchronous machines.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Terminals of a `Machine`
    enum class MachineTerminals : size_t
    {
      bus, ///< Bus of the machine, phases a, b and c
    };

    /// Inputs of a `Machine`, in per unit and of size 1
    enum class MachineInputs : size_t
    {
      efd, ///< \f$E_\mathrm{fd}\f$ Field voltage on the exciter base
      pm,  ///< \f$P_m\f$ Mechanical power
    };

    /// Parameters of a `Machine`
    enum class MachineParameters : size_t
    {
      S,    ///< \f$S_\mathrm{b}\f$ Rated three-phase apparent power [VA]
      V,    ///< \f$V_\mathrm{b}\f$ Rated line-to-line RMS voltage [V]
      f,    ///< \f$f_\mathrm{b}\f$ Rated electrical frequency [Hz]
      H,    ///< \f$H\f$ Inertia constant [s]
      F,    ///< \f$F\f$ Friction torque factor [p.u.]
      Rs,   ///< \f$R_s\f$ Stator winding resistance [p.u.]
      Ll,   ///< \f$L_l\f$ Stator leakage inductance [p.u.]
      Lmd,  ///< \f$L_\mathrm{md}\f$ Unsaturated d-axis magnetizing inductance [p.u.]
      Lmq,  ///< \f$L_\mathrm{mq}\f$ Unsaturated q-axis magnetizing inductance [p.u.]
      L0,   ///< \f$L_0\f$ Zero-sequence inductance [p.u.]
      Rfd,  ///< \f$R_\mathrm{fd}\f$ Field winding resistance [p.u.]
      Llfd, ///< \f$L_\mathrm{lfd}\f$ Field leakage inductance [p.u.]
      R1d,  ///< \f$R_\mathrm{1d}\f$ d-axis damper resistance [p.u.]
      Ll1d, ///< \f$L_\mathrm{l1d}\f$ d-axis damper leakage inductance [p.u.]
      R1q,  ///< \f$R_\mathrm{1q}\f$ q-axis damper 1 resistance [p.u.]
      Ll1q, ///< \f$L_\mathrm{l1q}\f$ q-axis damper 1 leakage inductance [p.u.]
      R2q,  ///< \f$R_\mathrm{2q}\f$ q-axis damper 2 resistance [p.u.]
      Ll2q, ///< \f$L_\mathrm{l2q}\f$ q-axis damper 2 leakage inductance [p.u.]
      S10,  ///< \f$S(1.0)\f$ Saturation factor at 1.0 per-unit flux
      S12,  ///< \f$S(1.2)\f$ Saturation factor at 1.2 per-unit flux
      Rfe,  ///< \f$R_\mathrm{fe}\f$ Iron-loss resistance [p.u.]
    };

    /// Data of a `Machine`
    template <typename RealT, typename IdxT>
    struct MachineData
    {
      using Parameters = MachineParameters;

      RealT S{ZERO<RealT>};
      RealT V{ZERO<RealT>};
      RealT f{ZERO<RealT>};
      RealT H{ZERO<RealT>};
      RealT F{ZERO<RealT>};
      RealT Rs{ZERO<RealT>};
      RealT Ll{ZERO<RealT>};
      RealT Lmd{ZERO<RealT>};
      RealT Lmq{ZERO<RealT>};
      RealT L0{ZERO<RealT>}; ///< Zero selects Ll
      RealT Rfd{ZERO<RealT>};
      RealT Llfd{ZERO<RealT>};
      RealT R1d{ZERO<RealT>};
      RealT Ll1d{ZERO<RealT>};
      RealT R1q{ZERO<RealT>};
      RealT Ll1q{ZERO<RealT>};
      RealT R2q{ZERO<RealT>};
      RealT Ll2q{ZERO<RealT>};
      RealT S10{ZERO<RealT>};
      RealT S12{ZERO<RealT>};
      RealT Rfe{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::S, &MachineData::S},
                                             Parameter{Parameters::V, &MachineData::V},
                                             Parameter{Parameters::f, &MachineData::f},
                                             Parameter{Parameters::H, &MachineData::H},
                                             Parameter{Parameters::F, &MachineData::F},
                                             Parameter{Parameters::Rs, &MachineData::Rs},
                                             Parameter{Parameters::Ll, &MachineData::Ll},
                                             Parameter{Parameters::Lmd, &MachineData::Lmd},
                                             Parameter{Parameters::Lmq, &MachineData::Lmq},
                                             Parameter{Parameters::L0, &MachineData::L0},
                                             Parameter{Parameters::Rfd, &MachineData::Rfd},
                                             Parameter{Parameters::Llfd, &MachineData::Llfd},
                                             Parameter{Parameters::R1d, &MachineData::R1d},
                                             Parameter{Parameters::Ll1d, &MachineData::Ll1d},
                                             Parameter{Parameters::R1q, &MachineData::R1q},
                                             Parameter{Parameters::Ll1q, &MachineData::Ll1q},
                                             Parameter{Parameters::R2q, &MachineData::R2q},
                                             Parameter{Parameters::Ll2q, &MachineData::Ll2q},
                                             Parameter{Parameters::S10, &MachineData::S10},
                                             Parameter{Parameters::S12, &MachineData::S12},
                                             Parameter{Parameters::Rfe, &MachineData::Rfe}};
    };
  } // namespace EMT
} // namespace GridKit

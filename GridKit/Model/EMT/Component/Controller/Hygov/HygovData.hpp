/**
 * @file HygovData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT hydro turbine governor.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `Hygov`, each of size 1 and in per unit
    enum class HygovInputs : size_t
    {
      omega, ///< \f$\omega\f$ Rotor speed
      pref,  ///< \f$P^\mathrm{ref}\f$ Load reference on the system power base
      paux,  ///< \f$P^\mathrm{aux}\f$ Auxiliary reference on the system power base
    };

    /// Parameters of a `Hygov`
    enum class HygovParameters : size_t
    {
      Rperm, ///< \f$R_{\mathrm{perm}}\f$ Permanent droop [p.u.]
      Rtemp, ///< \f$R_{\mathrm{temp}}\f$ Temporary droop [p.u.]
      Tr,    ///< \f$T_\mathrm{r}\f$ Temporary-droop reset time constant [s]
      Tf,    ///< \f$T_\mathrm{f}\f$ Governor error filter time constant [s]
      Tg,    ///< \f$T_\mathrm{g}\f$ Gate servo time constant [s]
      Velm,  ///< \f$V_{\mathrm{elm}}\f$ Maximum desired-gate velocity magnitude [p.u./s]
      Gmax,  ///< \f$G^{\max}\f$ Upper gate limit [p.u.]
      Gmin,  ///< \f$G^{\min}\f$ Lower gate limit [p.u.]
      Tw,    ///< \f$T_\mathrm{w}\f$ Water inertia time constant [s]
      At,    ///< \f$A_\mathrm{t}\f$ Turbine gain [p.u.]
      Dturb, ///< \f$D_{\mathrm{turb}}\f$ Turbine damping coefficient [p.u.]
      Qnl,   ///< \f$q_{\mathrm{NL}}\f$ No-load flow at nominal head [p.u.]
      Tn,    ///< \f$T_\mathrm{n}\f$ Speed lead-lag numerator time constant [s]
      Tnp,   ///< \f$T_{\mathrm{np}}\f$ Speed lead-lag denominator time constant [s]
      db1,   ///< \f$D_{\omega}\f$ Type 1 speed deadband threshold [p.u.]
      db2,   ///< \f$D_2\f$ Mechanical backlash, which is not supported [p.u.]
      Hdam,  ///< \f$H_{\mathrm{dam}}\f$ Dam head [p.u.]
      Gv0,   ///< \f$G_V^{(0)}\f$ Gate point 0 of the gain curve [p.u.]
      Gv1,   ///< \f$G_V^{(1)}\f$ Gate point 1 of the gain curve [p.u.]
      Gv2,   ///< \f$G_V^{(2)}\f$ Gate point 2 of the gain curve [p.u.]
      Gv3,   ///< \f$G_V^{(3)}\f$ Gate point 3 of the gain curve [p.u.]
      Gv4,   ///< \f$G_V^{(4)}\f$ Gate point 4 of the gain curve [p.u.]
      Gv5,   ///< \f$G_V^{(5)}\f$ Gate point 5 of the gain curve [p.u.]
      Pgv0,  ///< \f$P_{\mathrm{GV}}^{(0)}\f$ Power point 0 of the gain curve [p.u.]
      Pgv1,  ///< \f$P_{\mathrm{GV}}^{(1)}\f$ Power point 1 of the gain curve [p.u.]
      Pgv2,  ///< \f$P_{\mathrm{GV}}^{(2)}\f$ Power point 2 of the gain curve [p.u.]
      Pgv3,  ///< \f$P_{\mathrm{GV}}^{(3)}\f$ Power point 3 of the gain curve [p.u.]
      Pgv4,  ///< \f$P_{\mathrm{GV}}^{(4)}\f$ Power point 4 of the gain curve [p.u.]
      Pgv5,  ///< \f$P_{\mathrm{GV}}^{(5)}\f$ Power point 5 of the gain curve [p.u.]
      Trate, ///< \f$T^\mathrm{rate}\f$ Turbine rating [MW]
      S,     ///< \f$S\f$ Rated three-phase apparent power of the machine [VA]
    };

    /// Data of a `Hygov`
    template <typename RealT, typename IdxT>
    struct HygovData
    {
      using Parameters = HygovParameters;

      RealT Rperm{0.04};
      RealT Rtemp{0.3};
      RealT Tr{5.0};
      RealT Tf{0.05};
      RealT Tg{0.5};
      RealT Velm{0.2};
      RealT Gmax{ONE<RealT>};
      RealT Gmin{ZERO<RealT>};
      RealT Tw{ONE<RealT>};
      RealT At{1.2};
      RealT Dturb{0.5};
      RealT Qnl{0.05};
      RealT Tn{ZERO<RealT>};
      RealT Tnp{ZERO<RealT>};
      RealT db1{ZERO<RealT>};
      RealT db2{ZERO<RealT>};
      RealT Hdam{ONE<RealT>};
      RealT Gv0{ZERO<RealT>};
      RealT Gv1{ZERO<RealT>};
      RealT Gv2{ZERO<RealT>};
      RealT Gv3{ZERO<RealT>};
      RealT Gv4{ZERO<RealT>};
      RealT Gv5{ZERO<RealT>};
      RealT Pgv0{ZERO<RealT>};
      RealT Pgv1{ZERO<RealT>};
      RealT Pgv2{ZERO<RealT>};
      RealT Pgv3{ZERO<RealT>};
      RealT Pgv4{ZERO<RealT>};
      RealT Pgv5{ZERO<RealT>};
      RealT Trate{ZERO<RealT>};
      RealT S{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::Rperm, &HygovData::Rperm},
                                             Parameter{Parameters::Rtemp, &HygovData::Rtemp},
                                             Parameter{Parameters::Tr, &HygovData::Tr},
                                             Parameter{Parameters::Tf, &HygovData::Tf},
                                             Parameter{Parameters::Tg, &HygovData::Tg},
                                             Parameter{Parameters::Velm, &HygovData::Velm},
                                             Parameter{Parameters::Gmax, &HygovData::Gmax},
                                             Parameter{Parameters::Gmin, &HygovData::Gmin},
                                             Parameter{Parameters::Tw, &HygovData::Tw},
                                             Parameter{Parameters::At, &HygovData::At},
                                             Parameter{Parameters::Dturb, &HygovData::Dturb},
                                             Parameter{Parameters::Qnl, &HygovData::Qnl},
                                             Parameter{Parameters::Tn, &HygovData::Tn},
                                             Parameter{Parameters::Tnp, &HygovData::Tnp},
                                             Parameter{Parameters::db1, &HygovData::db1},
                                             Parameter{Parameters::db2, &HygovData::db2},
                                             Parameter{Parameters::Hdam, &HygovData::Hdam},
                                             Parameter{Parameters::Gv0, &HygovData::Gv0},
                                             Parameter{Parameters::Gv1, &HygovData::Gv1},
                                             Parameter{Parameters::Gv2, &HygovData::Gv2},
                                             Parameter{Parameters::Gv3, &HygovData::Gv3},
                                             Parameter{Parameters::Gv4, &HygovData::Gv4},
                                             Parameter{Parameters::Gv5, &HygovData::Gv5},
                                             Parameter{Parameters::Pgv0, &HygovData::Pgv0},
                                             Parameter{Parameters::Pgv1, &HygovData::Pgv1},
                                             Parameter{Parameters::Pgv2, &HygovData::Pgv2},
                                             Parameter{Parameters::Pgv3, &HygovData::Pgv3},
                                             Parameter{Parameters::Pgv4, &HygovData::Pgv4},
                                             Parameter{Parameters::Pgv5, &HygovData::Pgv5},
                                             Parameter{Parameters::Trate, &HygovData::Trate},
                                             Parameter{Parameters::S, &HygovData::S}};
    };
  } // namespace EMT
} // namespace GridKit

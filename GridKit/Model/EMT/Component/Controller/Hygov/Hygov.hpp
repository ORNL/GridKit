/**
 * @file Hygov.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT hydro turbine governor.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Hygov/HygovData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Hygov`, each of size 1 and in per unit
    enum class HygovInternalVariables : size_t
    {
      XN,      ///< \f$x_n\f$ Speed lead-lag denominator state
      XF,      ///< \f$x_f\f$ Governor error filter output on the turbine rating
      C,       ///< \f$c\f$ Desired-gate position on the turbine rating
      G,       ///< \f$g\f$ Gate position on the turbine rating
      Q,       ///< \f$q\f$ Turbine flow on the turbine rating
      OMEGADB, ///< \f$\omega_{\mathrm{db}}\f$ Deadbanded speed deviation
      EF,      ///< \f$e_f\f$ Governor error on the turbine rating
      FC,      ///< \f$f_c\f$ Desired-gate derivative target
      RC,      ///< \f$r_c\f$ Rate-limited desired-gate derivative target
      PGV,     ///< \f$P_{\mathrm{GV}}\f$ Gate-to-power curve output on the turbine rating
      H,       ///< \f$H\f$ Turbine head on the turbine rating
      PMECH,   ///< \f$P_{\mathrm{m}}\f$ Mechanical power on the machine rating
    };

    /// Hydro turbine governor
    template <typename scalar_type, typename index_type>
    class Hygov final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = HygovData<RealT, IdxT>;
      using InternalVariablesT = HygovInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = HygovInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Hygov";

      static constexpr std::array<IdxT, Utilities::enum_size<HygovInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>             NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<HygovInputs>()>            NUM_INPUT_VARIABLES{1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>               NUM_OUTPUT_VARIABLES{};

      explicit Hygov(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in HygovInternalVariables order
      enum : size_t
      {
        XN,
        XF,
        C,
        G,
        Q,
        OMEGADB,
        EF,
        FC,
        RC,
        PGV,
        H,
        PMECH,
      };

      // Offsets of the external values, in HygovInputs order
      enum : size_t
      {
        OMEGA,
        PREF,
        PAUX,
      };

      static constexpr RealT                           TIME_CONSTANT_MINIMUM{1.0e-3}; ///< Floor of `Tr`, `Tf`, `Tg`, `Tw` and `Tnp` [s]
      static constexpr size_t                          CURVE_POINTS{6};
      static constexpr std::array<RealT, CURVE_POINTS> DEFAULT_CURVE{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

      template <typename T>
      [[gnu::always_inline]] T gatePower(T gate) const;

      RealT Rperm_{ZERO<RealT>};
      RealT Rtemp_{ZERO<RealT>};
      RealT Tr_{ZERO<RealT>};
      RealT Tf_{ZERO<RealT>};
      RealT Tg_{ZERO<RealT>};
      RealT Velm_{ZERO<RealT>};
      RealT Gmax_{ZERO<RealT>};
      RealT Gmin_{ZERO<RealT>};
      RealT Tw_{ZERO<RealT>};
      RealT At_{ZERO<RealT>};
      RealT Dturb_{ZERO<RealT>};
      RealT Qnl_{ZERO<RealT>};
      RealT Tn_{ZERO<RealT>};
      RealT Tnp_{ZERO<RealT>};
      RealT db1_{ZERO<RealT>};
      RealT db2_{ZERO<RealT>};
      RealT Hdam_{ZERO<RealT>};
      RealT Trate_{ZERO<RealT>};
      RealT S_{ZERO<RealT>};

      std::array<RealT, CURVE_POINTS> Gv_{};  ///< Gate points of the gain curve
      std::array<RealT, CURVE_POINTS> Pgv_{}; ///< Power points of the gain curve

      RealT leadlag_gain_{ZERO<RealT>};       ///< Lead over lag time constant of the speed lead-lag
      RealT machine_to_turbine_{ZERO<RealT>}; ///< Machine rating over turbine rating
    };
  } // namespace EMT
} // namespace GridKit

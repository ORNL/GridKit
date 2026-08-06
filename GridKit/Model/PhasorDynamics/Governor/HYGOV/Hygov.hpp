/**
 * @file Hygov.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the HYGOV governor model.
 */

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class SignalNode;

    namespace Governor
    {
      /// Internal variables of a `Hygov`.
      enum class HygovInternalVariables : size_t
      {
        XN,      ///< \f$x_n\f$ Speed lead-lag denominator state
        XF,      ///< \f$x_f\f$ Governor error filter output on component base
        C,       ///< \f$c\f$ Desired-gate position on component base
        G,       ///< \f$g\f$ Gate position on component base
        Q,       ///< \f$q\f$ Turbine flow on component base
        OMEGADB, ///< \f$\omega_{\mathrm{db}}\f$ Deadbanded speed deviation
        EF,      ///< \f$e_f\f$ Governor error on component base
        FC,      ///< \f$f_c\f$ Desired-gate derivative target
        RC,      ///< \f$r_c\f$ Rate-limited desired-gate derivative target
        PGV,     ///< \f$P_{\mathrm{GV}}\f$ Gate-to-power curve output on component base
        H,       ///< \f$H\f$ Turbine head on component base
        PMECH,   ///< \f$P_{\mathrm{m}}\f$ Mechanical-power output on system base
        MAXIMUM,
      };

      /// External variables of a `Hygov`.
      enum class HygovExternalVariables : size_t
      {
        OMEGA, ///< \f$\omega\f$ Machine speed deviation
        PREF,  ///< \f$P^{\mathrm{ref}}\f$ Active-power/load reference on system base
        PAUX,  ///< \f$P^{\mathrm{aux}}\f$ Auxiliary power input on system base
        MAXIMUM,
      };

      /**
       * @brief Hydro turbine-governor model with temporary droop, gate servo,
       *        and a nonlinear single-penstock turbine.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class Hygov : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::va_system_base_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using SignalT            = SignalNode<ScalarT, IdxT>;
        using ModelDataT         = HygovData<RealT, IdxT>;
        using MonitorT           = Model::VariableMonitor<Hygov, HygovData>;
        using InternalVariablesT = HygovInternalVariables;
        using ExternalVariablesT = HygovExternalVariables;

        Hygov();
        explicit Hygov(const ModelDataT& data);
        ~Hygov();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                HygovInternalVariables,
                                HygovExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT* y,
            const ScalarT* yp,
            const ScalarT* wb,
            const ScalarT* ws,
            ScalarT*       f);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        /// Evaluate the nonlinear gate-to-power curve as a fixed sum of
        /// smooth linear segments.
        __attribute__((always_inline)) inline ScalarT gatePower(ScalarT gate) const;

        /// Steady component-base mechanical power at a gate and dam head.
        RealT initialMechanicalPower(RealT gate, RealT Hdam) const;

        /// Bisect a bracketed initialization residual to machine rounding.
        template <typename FuncT>
        static RealT bisectInitialRoot(RealT a,
                                       RealT b,
                                       RealT fa,
                                       RealT fb,
                                       FuncT residual);

        /// Solve the gate at the configured dam head.
        RealT solveInitialGate(RealT pmech) const;

        /// Solve the dam head that reproduces mechanical power at Gv5.
        RealT solveInitialDamHead(RealT pmech) const;

        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        /// Accepted seed distance beyond the achievable-power range edge.
        static constexpr RealT INITIALIZATION_TOLERANCE =
            static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

        RealT                Rperm_{static_cast<RealT>(0.04)};
        RealT                Rtemp_{static_cast<RealT>(0.3)};
        RealT                Tr_{static_cast<RealT>(5.0)};
        RealT                Tf_{static_cast<RealT>(0.05)};
        RealT                Tg_{static_cast<RealT>(0.5)};
        RealT                Velm_{static_cast<RealT>(0.2)};
        RealT                Gmax_{ONE<RealT>};
        RealT                Gmin_{ZERO<RealT>};
        RealT                Tw_{ONE<RealT>};
        RealT                At_{static_cast<RealT>(1.2)};
        RealT                Dturb_{static_cast<RealT>(0.5)};
        RealT                Qnl_{static_cast<RealT>(0.05)};
        RealT                Tn_{ZERO<RealT>};
        RealT                Tnp_{ZERO<RealT>};
        RealT                db1_{ZERO<RealT>};
        RealT                db2_{ZERO<RealT>};
        RealT                Hdam_{ONE<RealT>};
        std::array<RealT, 6> Gv_{};
        std::array<RealT, 6> Pgv_{};

        RealT va_component_base_{ZERO<RealT>};
        RealT leadlag_gain_{ZERO<RealT>};

        IdxT parameter_error_count_{0};

        RealT   Gmin_response_{Gmin_};
        RealT   Gmax_response_{Gmax_};
        RealT   Hdam_eff_{Hdam_};
        ScalarT pref_set_{0};
        ScalarT paux_set_{0};

        ComponentSignals<ScalarT, IdxT, HygovInternalVariables, HygovExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

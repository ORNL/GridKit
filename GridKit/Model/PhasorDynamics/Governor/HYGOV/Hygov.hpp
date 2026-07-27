/**
 * @file Hygov.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the HYGOV governor model.
 */

#pragma once

#include <array>
#include <cstddef>
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
        XN,      ///< Speed lead-lag denominator state
        XF,      ///< Governor error filter output
        C,       ///< Desired-gate position
        G,       ///< Gate position
        Q,       ///< Turbine flow
        OMEGADB, ///< Deadbanded speed deviation
        EF,      ///< Governor error into the filter
        FC,      ///< Desired-gate derivative target
        RC,      ///< Rate-limited desired-gate derivative target
        PGV,     ///< Nonlinear gate-to-power curve output
        H,       ///< Turbine head
        PMECH,   ///< Mechanical-power output
        MAXIMUM,
      };

      /// External variables of a `Hygov`.
      enum class HygovExternalVariables : size_t
      {
        OMEGA, ///< Machine speed deviation
        PREF,  ///< Active-power/load reference
        PAUX,  ///< Auxiliary power input
        MAXIMUM,
      };

      /// Indices into the HYGOV state, derivative, and residual vectors.
      struct HygovIdx
      {
        static constexpr size_t XN      = static_cast<size_t>(HygovInternalVariables::XN);
        static constexpr size_t XF      = static_cast<size_t>(HygovInternalVariables::XF);
        static constexpr size_t C       = static_cast<size_t>(HygovInternalVariables::C);
        static constexpr size_t G       = static_cast<size_t>(HygovInternalVariables::G);
        static constexpr size_t Q       = static_cast<size_t>(HygovInternalVariables::Q);
        static constexpr size_t OMEGADB = static_cast<size_t>(HygovInternalVariables::OMEGADB);
        static constexpr size_t EF      = static_cast<size_t>(HygovInternalVariables::EF);
        static constexpr size_t FC      = static_cast<size_t>(HygovInternalVariables::FC);
        static constexpr size_t RC      = static_cast<size_t>(HygovInternalVariables::RC);
        static constexpr size_t PGV     = static_cast<size_t>(HygovInternalVariables::PGV);
        static constexpr size_t H       = static_cast<size_t>(HygovInternalVariables::H);
        static constexpr size_t PMECH   = static_cast<size_t>(HygovInternalVariables::PMECH);
        static constexpr size_t MAXIMUM = static_cast<size_t>(HygovInternalVariables::MAXIMUM);
      };

      /// Indices into the HYGOV external-signal buffers.
      struct HygovExt
      {
        static constexpr size_t OMEGA   = static_cast<size_t>(HygovExternalVariables::OMEGA);
        static constexpr size_t PREF    = static_cast<size_t>(HygovExternalVariables::PREF);
        static constexpr size_t PAUX    = static_cast<size_t>(HygovExternalVariables::PAUX);
        static constexpr size_t MAXIMUM = static_cast<size_t>(HygovExternalVariables::MAXIMUM);
      };

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
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using ModelDataT = HygovData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Hygov, HygovData>;

        Hygov();
        explicit Hygov(const ModelDataT& data);
        ~Hygov();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT) override final;
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
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        /// Evaluate the nonlinear gate-to-power curve as a fixed sum of
        /// smooth linear segments.
        ScalarT gatePower(ScalarT gate) const;

        /// Analytic slope of the smooth gate-to-power curve, used to stamp
        /// the Jacobian entry the Enzyme auto-sparsity pass drops.
        RealT gatePowerDerivative(RealT gate) const;

        /// Solve the steady gate position that reproduces a seeded
        /// component-base mechanical power at an initial speed deviation.
        RealT solveInitialGate(RealT pmech, RealT omega) const;

        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        static constexpr RealT TIME_CONSTANT_MINIMUM    = static_cast<RealT>(1.0e-3);
        static constexpr RealT INITIALIZATION_TOLERANCE = static_cast<RealT>(1.0e-10);

        RealT                Trate_{ZERO<RealT>};
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

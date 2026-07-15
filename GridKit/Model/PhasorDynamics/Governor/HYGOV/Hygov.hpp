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

      template <typename scalar_type, typename index_type>
      class Hygov : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::allocated_;
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
        Hygov(const ModelDataT& data);
        ~Hygov() override;

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
        void initModelParams(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        ScalarT gatePower(ScalarT gate) const;
        RealT   invertGatePower(RealT pgv) const;
        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        RealT                Trate_{0};
        RealT                Rperm_{0};
        RealT                Rtemp_{0};
        RealT                Tr_{0};
        RealT                Tf_{0};
        RealT                Tg_{0};
        RealT                Velm_{0};
        RealT                Gmax_{0};
        RealT                Gmin_{0};
        RealT                Tw_{0};
        RealT                At_{0};
        RealT                Dturb_{0};
        RealT                Qnl_{0};
        RealT                Tn_{0};
        RealT                Tnp_{0};
        RealT                leadlag_gain_{0};
        RealT                db1_{0};
        RealT                db2_{0};
        RealT                Hdam_{1};
        std::array<RealT, 6> Gv_{};
        std::array<RealT, 6> Pgv_{};

        RealT va_component_base_{0};

        int parameter_error_count_{0};

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

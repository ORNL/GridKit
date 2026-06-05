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
    template <class ScalarT, typename IdxT>
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
        YOMEGA,  ///< Lead-lag-conditioned speed signal
        EF,      ///< Governor error into the filter
        FC,      ///< Desired-gate derivative target
        RC,      ///< Rate-limited desired-gate derivative target
        PGV,     ///< Nonlinear gate-to-power curve output
        PGVSAFE, ///< Safe gate-to-power value
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

      template <class ScalarT, typename IdxT>
      class Hygov : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::gridkit_component_id_;
        using Component<ScalarT, IdxT>::J_;
        using Component<ScalarT, IdxT>::J_cols_buffer_;
        using Component<ScalarT, IdxT>::J_rows_buffer_;
        using Component<ScalarT, IdxT>::J_vals_buffer_;
        using Component<ScalarT, IdxT>::residual_indices_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::variable_indices_;
        using Component<ScalarT, IdxT>::wb_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using model_data_type = HygovData<RealT, IdxT>;
        using MonitorT        = Model::VariableMonitor<Hygov, HygovData>;

        Hygov();
        Hygov(const model_data_type& data);
        ~Hygov();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
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
            ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        void initModelParams(const model_data_type& data);
        void initializeMonitor();

        ScalarT gatePower(ScalarT gate) const;
        RealT   invertGatePower(RealT pgv) const;
        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toPmechBase(ScalarT value) const;

        RealT                Trate_{0};
        RealT                pmech_mva_base_{0};
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

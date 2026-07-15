/**
 * @file Esdc1a.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the ESDC1A exciter model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <class ScalarT, typename IdxT>
    class SignalNode;

    namespace Exciter
    {
      /// Internal variables of an `Esdc1a`.
      enum class Esdc1aInternalVariables : size_t
      {
        EFDP, ///< Field-voltage state before optional speed multiplier
        VC,   ///< Sensed compensated voltage
        VR,   ///< Voltage-regulator output
        VF,   ///< Stabilizing feedback output
        XLL,  ///< Lead-lag state
        EV,   ///< Voltage-regulator input error
        VLL,  ///< Lead-lag block output
        VHV,  ///< High-value gate output
        SE,   ///< Saturation coefficient
        VFE,  ///< Exciter feedback signal
        EFD,  ///< Field-voltage output
        MAXIMUM,
      };

      /// External variables of an `Esdc1a`.
      enum class Esdc1aExternalVariables : size_t
      {
        OMEGA, ///< Machine speed deviation
        VREF,  ///< Voltage-control reference
        VS,    ///< Stabilizer input signal
        VUEL,  ///< Under-excitation limiter input
        MAXIMUM,
      };

      template <class ScalarT, typename IdxT>
      class Esdc1a : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::abs_tol_;
        using Component<ScalarT, IdxT>::allocated_;
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::gridkit_component_id_;
        using Component<ScalarT, IdxT>::J_cols_buffer_;
        using Component<ScalarT, IdxT>::J_rows_buffer_;
        using Component<ScalarT, IdxT>::J_vals_buffer_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::residual_indices_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::variable_indices_;
        using Component<ScalarT, IdxT>::wb_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using bus_type        = BusBase<ScalarT, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using model_data_type = Esdc1aData<RealT, IdxT>;
        using MonitorT        = Model::VariableMonitor<Esdc1a, Esdc1aData>;

        Esdc1a(bus_type* bus);
        Esdc1a(bus_type* bus, const model_data_type& data);
        ~Esdc1a();

        int setGridKitComponentID(IdxT) override final;
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
                                Esdc1aInternalVariables,
                                Esdc1aExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        void initModelParams(const model_data_type& data);
        void setDerivedParams();
        void initializeMonitor();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        bus_type* bus_{nullptr};

        RealT Tr_{0.0};
        RealT Ka_{40.0};
        RealT Ta_{0.1};
        RealT Tb_{0.0};
        RealT Tc_{0.0};
        RealT Vrmax_{1.0};
        RealT Vrmin_{-1.0};
        RealT Ke_{0.1};
        RealT Te_{0.5};
        RealT Kf_{0.05};
        RealT Tf1_{0.7};
        RealT spdmlt_{0.0};
        RealT E1_{2.8};
        RealT Se1_{0.08};
        RealT E2_{3.7};
        RealT Se2_{0.33};
        IdxT  UEL_{0};
        RealT exclim_{1.0};

        IdxT parameter_error_count_{0};

        RealT sUEL_{0};
        RealT sUELoff_{1};
        RealT slim_{0};
        RealT slim_off_{1};
        RealT SA_{0};
        RealT SB_{0};

        ScalarT vref_{0};

        ComponentSignals<ScalarT, IdxT, Esdc1aInternalVariables, Esdc1aExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                         monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

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
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
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

      /// Indices into the ESDC1A state, derivative, and residual vectors.
      struct Esdc1aIdx
      {
        static constexpr size_t EFDP    = static_cast<size_t>(Esdc1aInternalVariables::EFDP);
        static constexpr size_t VC      = static_cast<size_t>(Esdc1aInternalVariables::VC);
        static constexpr size_t VR      = static_cast<size_t>(Esdc1aInternalVariables::VR);
        static constexpr size_t VF      = static_cast<size_t>(Esdc1aInternalVariables::VF);
        static constexpr size_t XLL     = static_cast<size_t>(Esdc1aInternalVariables::XLL);
        static constexpr size_t EV      = static_cast<size_t>(Esdc1aInternalVariables::EV);
        static constexpr size_t VLL     = static_cast<size_t>(Esdc1aInternalVariables::VLL);
        static constexpr size_t VHV     = static_cast<size_t>(Esdc1aInternalVariables::VHV);
        static constexpr size_t SE      = static_cast<size_t>(Esdc1aInternalVariables::SE);
        static constexpr size_t VFE     = static_cast<size_t>(Esdc1aInternalVariables::VFE);
        static constexpr size_t EFD     = static_cast<size_t>(Esdc1aInternalVariables::EFD);
        static constexpr size_t MAXIMUM = static_cast<size_t>(Esdc1aInternalVariables::MAXIMUM);
      };

      /// Indices into the ESDC1A external-signal buffers.
      struct Esdc1aExt
      {
        static constexpr size_t OMEGA   = static_cast<size_t>(Esdc1aExternalVariables::OMEGA);
        static constexpr size_t VREF    = static_cast<size_t>(Esdc1aExternalVariables::VREF);
        static constexpr size_t VS      = static_cast<size_t>(Esdc1aExternalVariables::VS);
        static constexpr size_t VUEL    = static_cast<size_t>(Esdc1aExternalVariables::VUEL);
        static constexpr size_t MAXIMUM = static_cast<size_t>(Esdc1aExternalVariables::MAXIMUM);
      };

      template <typename scalar_type, typename index_type>
      class Esdc1a : public Component<scalar_type, index_type>
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
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using BusT       = BusBase<ScalarT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using ModelDataT = Esdc1aData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Esdc1a, Esdc1aData>;

        Esdc1a(BusT* bus);
        Esdc1a(BusT* bus, const ModelDataT& data);
        ~Esdc1a();

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
                                Esdc1aInternalVariables,
                                Esdc1aExternalVariables>&
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

        /// Recover the input that the smooth CommonMath ramp maps to a
        /// requested strictly positive output.
        RealT inverseRamp(RealT ramp_output) const;

        ScalarT& Vr();
        ScalarT& Vi();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        BusT* bus_{nullptr};

        RealT Tr_{ZERO<RealT>};
        RealT Ka_{static_cast<RealT>(40.0)};
        RealT Ta_{static_cast<RealT>(0.1)};
        RealT Tb_{ZERO<RealT>};
        RealT Tc_{ZERO<RealT>};
        RealT Vrmax_{ONE<RealT>};
        RealT Vrmin_{static_cast<RealT>(-1.0)};
        RealT Ke_{static_cast<RealT>(0.1)};
        RealT Te_{static_cast<RealT>(0.5)};
        RealT Kf_{static_cast<RealT>(0.05)};
        RealT Tf1_{static_cast<RealT>(0.7)};
        bool  Spdmlt_{false};
        RealT E1_{static_cast<RealT>(2.8)};
        RealT Se1_{static_cast<RealT>(0.08)};
        RealT E2_{static_cast<RealT>(3.7)};
        RealT Se2_{static_cast<RealT>(0.33)};
        IdxT  UEL_{0};
        bool  exclim_{true};
        RealT spd_on_{0};
        RealT uel_on_{0};
        RealT uel_off_{1};
        RealT lim_on_{1};
        RealT lim_off_{0};
        RealT SA_{0};
        RealT SB_{0};

        IdxT parameter_error_count_{0};

        ScalarT omega_set_{0};
        ScalarT vref_set_{0};
        ScalarT vs_set_{0};
        ScalarT vuel_set_{0};

        ComponentSignals<ScalarT, IdxT, Esdc1aInternalVariables, Esdc1aExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                         monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

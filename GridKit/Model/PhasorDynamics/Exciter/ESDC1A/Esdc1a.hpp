/**
 * @file Esdc1a.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the ESDC1A exciter model.
 */

#pragma once

#include <cstddef>
#include <memory>

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
        EFDP,    ///< \f$E_{\mathrm{fd}}'\f$ Differential exciter field-voltage state [p.u.]
        VC,      ///< \f$V_C\f$ Differential filtered terminal-voltage magnitude [p.u.]
        VR,      ///< \f$V_R\f$ Differential voltage-regulator output [p.u.]
        VF,      ///< \f$V_F\f$ Differential stabilizing feedback state [p.u.]
        XLL,     ///< \f$x_{\mathrm{LL}}\f$ Differential input lead-lag denominator state [p.u.]
        EV,      ///< \f$e_V\f$ Algebraic voltage-error summing output [p.u.]
        VLL,     ///< \f$V_{\mathrm{LL}}\f$ Algebraic input lead-lag output [p.u.]
        VHV,     ///< \f$V_{\mathrm{HV}}\f$ Algebraic high-value gate output [p.u.]
        SE,      ///< \f$s_e\f$ Scaled-quadratic saturation contribution [p.u.]
        VFE,     ///< \f$V_{\mathrm{FE}}\f$ Algebraic exciter feedback drive [p.u.]
        EFD,     ///< \f$E_{\mathrm{fd}}\f$ Algebraic field-voltage output [p.u.]
        MAXIMUM, ///< Number of ESDC1A internal variables
      };

      /// External signal variables read or initialized by an `Esdc1a`.
      enum class Esdc1aExternalVariables : size_t
      {
        OMEGA,   ///< \f$\omega\f$ Known machine speed deviation [p.u.]
        VREF,    ///< \f$V_{\mathrm{ref}}\f$ Unknown voltage-control reference [p.u.]
        VS,      ///< \f$V_S\f$ Known stabilizer input signal [p.u.]
        VUEL,    ///< \f$V_{\mathrm{UEL}}\f$ Known under-excitation limiter input [p.u.]
        MAXIMUM, ///< Number of ESDC1A external signal variables
      };

      /**
       * @brief IEEE DC1A excitation-system model (ESDC1A).
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class Esdc1a : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::alpha;
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
        using Component<scalar_type, index_type>::ws_;
        using Component<scalar_type, index_type>::ws_indices_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using BusT               = BusBase<ScalarT, IdxT>;
        using SignalT            = SignalNode<ScalarT, IdxT>;
        using ModelDataT         = Esdc1aData<RealT, IdxT>;
        using MonitorT           = Model::VariableMonitor<Esdc1a, Esdc1aData>;
        using InternalVariablesT = Esdc1aInternalVariables;
        using ExternalVariablesT = Esdc1aExternalVariables;

        Esdc1a(BusT* bus);
        Esdc1a(BusT* bus, const ModelDataT& data);
        ~Esdc1a();

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
                                Esdc1aInternalVariables,
                                Esdc1aExternalVariables>&
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

        static __attribute__((always_inline)) inline ScalarT awmin(
            ScalarT x,
            ScalarT f,
            RealT   xmin);

        /// Recover the input that the smooth CommonMath ramp maps to a
        /// requested strictly positive output.
        RealT inverseRamp(RealT ramp_output) const;

        ScalarT& Vr();
        ScalarT& Vi();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static void            logTimeConstantWarning();

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
        RealT lim_on_{1};
        RealT SA_{0};
        RealT SB_{0};
        RealT Ke_eff_{Ke_};

        IdxT parameter_error_count_{0};

        ScalarT omega_set_{0};
        ScalarT vref_set_{0};
        ScalarT vs_set_{0};
        ScalarT vuel_set_{0};

        ComponentSignals<ScalarT, IdxT, Esdc1aInternalVariables, Esdc1aExternalVariables>
                                  signals_;
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

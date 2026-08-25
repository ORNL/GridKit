/**
 * @file Regca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REGCA phasor-dynamics converter model.
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <memory>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
    class SignalNode;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Internal variables of a `Regca`
      enum class RegcaInternalVariables : size_t
      {
        VM,      ///< \f$V_M\f$ Filtered terminal voltage
        IQ,      ///< \f$I_q\f$ Reactive-current state on component base
        IP,      ///< \f$I_p\f$ Active-current state on component base
        VT,      ///< \f$V_T\f$ Terminal voltage magnitude
        IR,      ///< \f$I_\mathrm{r}\f$ Branch-current real component on system base
        II,      ///< \f$I_\mathrm{i}\f$ Branch-current imaginary component on system base
        IQEXTRA, ///< \f$I_q^\mathrm{extra}\f$ Extra inductive current from HVRCM on component base
        IL,      ///< \f$I_L\f$ LVPL upper-limit current curve on component base
        PBR,     ///< \f$P^\mathrm{br}\f$ Branch active power on system base
        QBR,     ///< \f$Q^\mathrm{br}\f$ Branch reactive power on system base
        MAXIMUM,
      };

      /// External variables of a `Regca`
      enum class RegcaExternalVariables : size_t
      {
        IPCMD, ///< \f$I_p^\mathrm{cmd}\f$ Active-current command on system base
        IQCMD, ///< \f$I_q^\mathrm{cmd}\f$ Reactive-current command on system base
        MAXIMUM,
      };

      /**
       * @brief First-generation WECC renewable generator/converter model (REGCA).
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class Regca : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::h_;
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
        using BusT               = BusBase<ScalarT, IdxT>;
        using SignalT            = SignalNode<ScalarT, IdxT>;
        using ModelDataT         = RegcaData<RealT, IdxT>;
        using MonitorT           = Model::VariableMonitor<Regca, RegcaData>;
        using InternalVariablesT = RegcaInternalVariables;
        using ExternalVariablesT = RegcaExternalVariables;

        Regca(BusT* bus);
        Regca(BusT* bus, const ModelDataT& data);
        ~Regca();

        int  setGridKitComponentID(IdxT component_id) override final;
        int  allocate() override final;
        int  verify() const override final;
        int  initialize() override final;
        int  tagDifferentiable() override final;
        int  setAbsoluteTolerance(RealT rel_tol) override final;
        void setLvplGain(RealT KL);
        int  evaluateResidual() override final;
        int  evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                RegcaInternalVariables,
                                RegcaExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* ws, ScalarT* f);

        __attribute__((always_inline)) inline int evaluateBusResidual(
            const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* h);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        /**
         * @brief Smooth approximation of the REGCA `rrpwr` rate limiter.
         *
         * Limits motion that increases \f$|x|\f$ while smoothly releasing
         * restoring motion toward zero. At \f$x = 0\f$, the result is the
         * symmetric smooth slew limit applied to @p f.
         *
         * @param[in] x State whose sign determines the limited direction.
         * @param[in] f Unconstrained derivative.
         * @param[in] rate Nonnegative symmetric rate limit.
         * @return Rate-limited derivative.
         *
         * @todo Move this reusable limiter to CommonMath.
         */
        static __attribute__((always_inline)) inline ScalarT rrpwr(
            const ScalarT x,
            const ScalarT f,
            const RealT   rate)
        {
          assert(rate >= ZERO<RealT>);

          const ScalarT t     = std::tanh(HALF<RealT> * Math::MU<RealT> * x);
          const ScalarT abs_t = std::abs(t);
          const ScalarT w_pos = HALF<RealT> * t * (t + abs_t);
          const ScalarT w_neg = HALF<RealT> * t * (t - abs_t);

          return f + (ONE<RealT> - w_pos) * Math::ramp(-f - rate)
                 - (ONE<RealT> - w_neg) * Math::ramp(f - rate);
        }

        /**
         * @brief Smooth anti-windup derivative under a moving upper bound.
         *
         * Below the bound the unconstrained derivative passes. Pinned at or
         * above the bound, the state tracks min(f, fmax), so a falling bound
         * drags the state down with it. Equivalent to fixed-bound anti-windup
         * on the gap x - xmax held below zero.
         *
         * @param[in] x State limited from above.
         * @param[in] f Unconstrained derivative of x.
         * @param[in] xmax Moving upper bound on x.
         * @param[in] fmax Derivative of the moving upper bound.
         * @return Anti-windup-limited derivative.
         *
         * @todo Move this one-sided anti-windup helper to CommonMath.
         */
        static __attribute__((always_inline)) inline ScalarT awmax(
            const ScalarT x,
            const ScalarT f,
            const ScalarT xmax,
            const ScalarT fmax)
        {
          const ScalarT gap_rate = f - fmax;
          const ScalarT below    = Math::sigmoid(xmax - x);

          return fmax
                 + (below + (ONE<RealT> - below) * Math::sigmoid(-gap_rate))
                       * gap_rate;
        }

        ScalarT initialHvrcmCurrent(ScalarT dv) const;

        ScalarT& Vr();
        ScalarT& Vi();
        ScalarT& Ir();
        ScalarT& Ii();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static void logTimeConstantWarning();

        BusT* bus_{nullptr};

        // Input parameters
        RealT p0_{0};
        RealT q0_{0};
        RealT mva_base_{0};
        RealT Tg_{0};
        RealT TM_{0};
        RealT Rqmax_{0};
        RealT Rqmin_{0};
        RealT Rpmax_{0};
        bool  sL_{false};
        RealT IL1_{0};
        RealT VL0_{0};
        RealT VL1_{0};
        RealT VA0_{0};
        RealT VA1_{0};
        RealT Vhvmax_{0};
        RealT KL_{100.0};

        IdxT parameter_error_count_{0};

        // Derived parameters
        RealT va_converter_base_{0};
        RealT use_lvpl_{0};
        RealT bypass_lvpl_{1};
        RealT use_rqmax_{0};
        RealT use_rqmin_{0};

        // Unattached command setpoints
        ScalarT ipcmd_set_{0};
        ScalarT iqcmd_set_{0};

        ComponentSignals<ScalarT, IdxT, RegcaInternalVariables, RegcaExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        // Local copies of signal variables
        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

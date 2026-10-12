/**
 * @file Ieeest.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the IEEEST Power System Stabilizer.
 */

#pragma once

#include <cstddef>
#include <memory>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Internal variables of an `Ieeest`. Notch states come last; a model
      /// of notch order \f$n\f$ has the first \f$n\f$ of them.
      enum class IeeestInternalVariables : size_t
      {
        X5,  ///< \f$x_5\f$ Lead-lag 1 state
        X6,  ///< \f$x_6\f$ Lead-lag 2 state
        X7,  ///< \f$x_7\f$ Washout state
        VSS, ///< \f$V_{\mathrm{ss}}\f$ Limited stabilizer signal and model output
        X1,  ///< \f$x_1\f$ First notch factor output
        X2,  ///< \f$x_2\f$ Rate of \f$x_1\f$
        X3,  ///< \f$x_3\f$ Second notch factor output
        X4,  ///< \f$x_4\f$ Rate of \f$x_3\f$
      };

      /// External variables of an `Ieeest`.
      enum class IeeestExternalVariables : size_t
      {
        U, ///< \f$u\f$ Stabilizer input signal
      };

      /// Degree of \f$1 + p s + q s^2\f$. A zero coefficient marks an absent term.
      template <typename RealT>
      constexpr size_t factorDegree(RealT p, RealT q)
      {
        if (q != RealT{0})
        {
          return 2;
        }
        if (p != RealT{0})
        {
          return 1;
        }
        return 0;
      }

      /// Degree of the notch denominator
      /// \f$(1 + A_1 s + A_2 s^2)(1 + A_3 s + A_4 s^2)\f$.
      template <typename RealT>
      constexpr size_t notchOrder(RealT A1, RealT A2, RealT A3, RealT A4)
      {
        return factorDegree(A1, A2) + factorDegree(A3, A4);
      }

      /**
       * @brief IEEE power system stabilizer (IEEEST).
       *
       * IeeestFactory selects the notch order from the model data.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       * @tparam order Degree of the notch denominator, at most four.
       */
      template <typename scalar_type, typename index_type, size_t order>
      class Ieeest : public Component<scalar_type, index_type>
      {
        static_assert(order <= 4, "IEEEST notch order is at most four");

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
        using Component<scalar_type, index_type>::ws_;
        using Component<scalar_type, index_type>::ws_indices_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT         = IeeestData<RealT, IdxT>;
        using SignalNodeSetT     = SignalNodeSet<ScalarT, IdxT>;
        using SignalPortsT       = SignalPorts<ScalarT, ModelDataT>;
        using MonitorT           = Model::VariableMonitor<Ieeest, IeeestData>;
        using InternalVariablesT = IeeestInternalVariables;
        using ExternalVariablesT = IeeestExternalVariables;

        Ieeest();
        Ieeest(const ModelDataT& data);
        ~Ieeest();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        Model::ConfigurationChecks verify() const override final;

        SignalPortsT& getPorts()
        {
          return ports_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        [[gnu::always_inline]] inline int evaluateInternalResidual(
            const ScalarT*,
            const ScalarT*,
            const ScalarT*,
            const ScalarT*,
            ScalarT*);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static void            logTimeConstantWarning();

        RealT A1_{0};
        RealT A2_{0};
        RealT A3_{0};
        RealT A4_{0};
        RealT A5_{0};
        RealT A6_{0};
        RealT T1_{0};
        RealT T2_{1};
        RealT T3_{0};
        RealT T4_{1};
        RealT T5_{0};
        RealT T6_{1};
        RealT Ks_{1};
        RealT Lsmin_{-0.1};
        RealT Lsmax_{0.1};
        RealT Vcl_{0};
        RealT Vcu_{0};
        RealT Tdelay_{0};

        RealT b1_{0}; ///< First notch factor \f$1 + b_1 s + b_2 s^2\f$, which carries the numerator
        RealT b2_{0};
        RealT c1_{0}; ///< Second notch factor \f$1 + c_1 s + c_2 s^2\f$
        RealT c2_{0};

        SignalPortsT              ports_;
        std::unique_ptr<MonitorT> monitor_;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit

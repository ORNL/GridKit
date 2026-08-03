/**
 * @file Ieeest.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the IEEEST Power System Stabilizer.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      template <typename real_type, typename index_type>
      struct IeeestData;
    } // namespace Stabilizer

    template <typename scalar_type, typename index_type>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Internal variables of a `Ieeest`
      enum class IeeestInternalVariables : size_t
      {
        X1,  ///< \f$x_1\f$ Notch filter signal state [p.u.], differential for \f$n\ge1\f$
        X2,  ///< \f$x_2\f$ First derivative of the filtered signal [p.u./sec], differential for \f$n\ge2\f$
        X3,  ///< \f$x_3\f$ Second derivative of the filtered signal [p.u./sec^2], differential for \f$n\ge3\f$
        X4,  ///< \f$x_4\f$ Third derivative of the filtered signal [p.u./sec^3], differential for \f$n=4\f$
        X5,  ///< \f$x_5\f$ Lead-lag 1 state [p.u.]
        X6,  ///< \f$x_6\f$ Lead-lag 2 state [p.u.]
        X7,  ///< \f$x_7\f$ Washout state [p.u.]
        V4,  ///< \f$v_4\f$ Notch filter output [p.u.]
        V5,  ///< \f$v_5\f$ Lead-lag 1 output [p.u.]
        V6,  ///< \f$v_6\f$ Lead-lag 2 output [p.u.]
        V7,  ///< \f$v_7\f$ Unlimited stabilizer signal [p.u.]
        VSS, ///< \f$V_{ss}\f$ Limited stabilizer signal, the model output [p.u.]
        MAXIMUM,
      };

      /// External variables of a `Ieeest`
      enum class IeeestExternalVariables : size_t
      {
        U, ///< \f$u\f$ Stabilizer input signal [p.u.]
        MAXIMUM,
      };

      /**
       * @brief IEEE type ST power system stabilizer (IEEEST).
       *
       * A selectable-order notch filter, two lead-lag blocks, a washout, and an
       * output limiter. The notch order \f$n\in\{0,1,2,3,4\}\f$ is the degree of
       * the expanded denominator derived from \f$A_1,\ldots,A_4\f$.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       *
       * @see IeeestData
       */
      template <typename scalar_type, typename index_type>
      class Ieeest : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::time_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::h_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::allocated_;

      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT = IeeestData<RealT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Ieeest, IeeestData>;

        Ieeest();
        Ieeest(const ModelDataT& data);
        ~Ieeest();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        /// Get the `ComponentSignals` from this `Ieeest`
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                IeeestInternalVariables,
                                IeeestExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*,
            const ScalarT*,
            const ScalarT*,
            const ScalarT*,
            ScalarT*);

      private:
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

        /// Expanded notch denominator \f$a_0 + a_1 s + a_2 s^2 + a_3 s^3 + a_4 s^4\f$
        RealT a0_{1};
        RealT a1_{0};
        RealT a2_{0};
        RealT a3_{0};
        RealT a4_{0};

        /// Order indicators \f$O_k\f$, one when the notch order is at least \f$k\f$
        RealT O1_{0};
        RealT O2_{0};
        RealT O3_{0};

        // A denominator time constant that is not positive bypasses its block.
        RealT bypass_T2_{0};
        RealT bypass_T4_{0};
        RealT bypass_T6_{0};

        ComponentSignals<ScalarT, IdxT, IeeestInternalVariables, IeeestExternalVariables> signals_;

        std::unique_ptr<MonitorT> monitor_;

        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit

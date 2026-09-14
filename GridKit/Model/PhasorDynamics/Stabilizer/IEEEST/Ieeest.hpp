/**
 * @file Ieeest.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the IEEEST Power System Stabilizer.
 */

#pragma once

#include <cstddef>
#include <memory>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class SignalNode;

    namespace Stabilizer
    {
      /// Complete internal variable layout for each notch-filter order.
      template <size_t order>
      struct IeeestVariables;

      template <>
      struct IeeestVariables<0>
      {
        enum class InternalVariables : size_t
        {
          X5,  ///< Lead-lag 1 state
          X6,  ///< Lead-lag 2 state
          X7,  ///< Washout state
          V4,  ///< Notch-filter output
          V5,  ///< Lead-lag 1 output
          V6,  ///< Lead-lag 2 output
          V7,  ///< Unlimited stabilizer signal
          VSS, ///< Limited stabilizer signal and model output
        };
      };

      template <>
      struct IeeestVariables<1>
      {
        enum class InternalVariables : size_t
        {
          X1,  ///< Notch-filter state
          X5,  ///< Lead-lag 1 state
          X6,  ///< Lead-lag 2 state
          X7,  ///< Washout state
          V4,  ///< Notch-filter output
          V5,  ///< Lead-lag 1 output
          V6,  ///< Lead-lag 2 output
          V7,  ///< Unlimited stabilizer signal
          VSS, ///< Limited stabilizer signal and model output
        };
      };

      template <>
      struct IeeestVariables<2>
      {
        enum class InternalVariables : size_t
        {
          X1,  ///< Notch-filter state
          X2,  ///< First derivative of X1
          X5,  ///< Lead-lag 1 state
          X6,  ///< Lead-lag 2 state
          X7,  ///< Washout state
          V4,  ///< Notch-filter output
          V5,  ///< Lead-lag 1 output
          V6,  ///< Lead-lag 2 output
          V7,  ///< Unlimited stabilizer signal
          VSS, ///< Limited stabilizer signal and model output
        };
      };

      template <>
      struct IeeestVariables<3>
      {
        enum class InternalVariables : size_t
        {
          X1,  ///< Notch-filter state
          X2,  ///< First derivative of X1
          X3,  ///< Second derivative of X1
          X5,  ///< Lead-lag 1 state
          X6,  ///< Lead-lag 2 state
          X7,  ///< Washout state
          V4,  ///< Notch-filter output
          V5,  ///< Lead-lag 1 output
          V6,  ///< Lead-lag 2 output
          V7,  ///< Unlimited stabilizer signal
          VSS, ///< Limited stabilizer signal and model output
        };
      };

      template <>
      struct IeeestVariables<4>
      {
        enum class InternalVariables : size_t
        {
          X1,  ///< Notch-filter state
          X2,  ///< First derivative of X1
          X3,  ///< Second derivative of X1
          X4,  ///< Third derivative of X1
          X5,  ///< Lead-lag 1 state
          X6,  ///< Lead-lag 2 state
          X7,  ///< Washout state
          V4,  ///< Notch-filter output
          V5,  ///< Lead-lag 1 output
          V6,  ///< Lead-lag 2 output
          V7,  ///< Unlimited stabilizer signal
          VSS, ///< Limited stabilizer signal and model output
        };
      };

      /// Internal variables of a `Ieeest` of the given notch-filter order
      template <size_t order>
      using IeeestInternalVariables = typename IeeestVariables<order>::InternalVariables;

      /// External variables of a `Ieeest`
      enum class IeeestExternalVariables : size_t
      {
        U, ///< \f$u\f$ Stabilizer input signal
      };

      /**
       * @brief IEEEST with a compile-time notch-filter order in [0, 4].
       *
       * Only the active notch states are allocated. IeeestFactory selects the
       * specialization from the two notch-denominator factors.
       *
       * @tparam scalar_type Scalar data type.
       * @tparam index_type Index data type.
       * @tparam order Degree of the expanded notch denominator.
       */
      template <typename scalar_type, typename index_type, size_t order>
      class Ieeest : public Component<scalar_type, index_type>
      {
        static_assert(order <= 4, "Ieeest notch filter order must be in [0, 4]");

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
        using Component<scalar_type, index_type>::ws_;
        using Component<scalar_type, index_type>::ws_indices_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using SignalT            = SignalNode<ScalarT, IdxT>;
        using ModelDataT         = IeeestData<RealT, IdxT>;
        using SignalNodeSetT     = SignalNodeSet<ScalarT, IdxT>;
        using SignalPortsT       = SignalPorts<ScalarT, ModelDataT>;
        using MonitorT           = Model::VariableMonitor<Ieeest, IeeestData>;
        using InternalVariablesT = IeeestInternalVariables<order>;
        using ExternalVariablesT = IeeestExternalVariables;

        Ieeest();
        explicit Ieeest(const ModelDataT& data);
        ~Ieeest();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

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
        static constexpr size_t X5  = static_cast<size_t>(InternalVariablesT::X5);
        static constexpr size_t X6  = static_cast<size_t>(InternalVariablesT::X6);
        static constexpr size_t X7  = static_cast<size_t>(InternalVariablesT::X7);
        static constexpr size_t V4  = static_cast<size_t>(InternalVariablesT::V4);
        static constexpr size_t V5  = static_cast<size_t>(InternalVariablesT::V5);
        static constexpr size_t V6  = static_cast<size_t>(InternalVariablesT::V6);
        static constexpr size_t V7  = static_cast<size_t>(InternalVariablesT::V7);
        static constexpr size_t VSS = static_cast<size_t>(InternalVariablesT::VSS);
        static constexpr size_t U   = static_cast<size_t>(IeeestExternalVariables::U);

        void loadRealParameter(const ModelDataT& data, IeeestParameters parameter, RealT& value);
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

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

        RealT a1_{0};
        RealT a2_{0};
        RealT a3_{0};
        RealT a4_{0};
        RealT inv_an_{0}; ///< Reciprocal of the active leading denominator coefficient.
        RealT inv_T2_{1};
        RealT inv_T4_{1};
        RealT inv_T6_{1};

        bool parameters_valid_{true};

        SignalPortsT              ports_;
        std::unique_ptr<MonitorT> monitor_;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit

/**
 * @file Ieeest.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the IEEEST Power System Stabilizer.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
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
      /// Internal variables of a `Ieeest`.
      enum class IeeestInternalVariables : size_t
      {
        X1,  ///< \f$x_1\f$ Notch-filter signal state
        X2,  ///< \f$x_2\f$ First derivative of the filtered signal
        X3,  ///< \f$x_3\f$ Second derivative of the filtered signal
        X4,  ///< \f$x_4\f$ Third derivative of the filtered signal
        X5,  ///< \f$x_5\f$ Lead-lag 1 state
        X6,  ///< \f$x_6\f$ Lead-lag 2 state
        X7,  ///< \f$x_7\f$ Washout state
        V4,  ///< \f$v_4\f$ Notch-filter output
        V5,  ///< \f$v_5\f$ Lead-lag 1 output
        V6,  ///< \f$v_6\f$ Lead-lag 2 output
        V7,  ///< \f$v_7\f$ Unlimited stabilizer signal
        VSS, ///< \f$V_{\mathrm{ss}}\f$ Limited stabilizer signal and model output
        MAXIMUM,
      };

      /// External variables of a `Ieeest`.
      enum class IeeestExternalVariables : size_t
      {
        U, ///< \f$u\f$ Stabilizer input signal
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class Ieeest : public Component<scalar_type, index_type>
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
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using SignalT            = SignalNode<ScalarT, IdxT>;
        using ModelDataT         = IeeestData<RealT, IdxT>;
        using MonitorT           = Model::VariableMonitor<Ieeest, IeeestData>;
        using InternalVariablesT = IeeestInternalVariables;
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
            const ScalarT* y,
            const ScalarT* yp,
            const ScalarT* wb,
            const ScalarT* ws,
            ScalarT*       f);

      private:
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

        IdxT order_{0};

        RealT a1_{0};
        RealT a2_{0};
        RealT a3_{0};
        RealT a4_{0};

        IdxT parameter_error_count_{0};

        ComponentSignals<ScalarT, IdxT, IeeestInternalVariables, IeeestExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                         monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit

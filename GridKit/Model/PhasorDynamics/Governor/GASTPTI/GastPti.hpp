/**
 * @file GastPti.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the GASTPTI governor model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class SignalNode;

    namespace Governor
    {
      /// Internal variables of a `GastPti`.
      enum class GastPtiInternalVariables : size_t
      {
        XVALVE, ///< Fuel-valve state
        XFLOW,  ///< Fuel-flow state
        XTEMP,  ///< Exhaust-temperature feedback state
        VLOAD,  ///< Speed/load fuel demand
        VTEMP,  ///< Temperature-limit fuel demand
        VLV,    ///< LV gate output
        PMECH,  ///< Mechanical power output
        MAXIMUM,
      };

      /// External variables of a `GastPti`.
      enum class GastPtiExternalVariables : size_t
      {
        OMEGA, ///< Machine speed deviation
        PREF,  ///< Active-power/load reference
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class GastPti : public Component<scalar_type, index_type>
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
        using Component<scalar_type, index_type>::va_system_base_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT = GastPtiData<RealT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using MonitorT   = Model::VariableMonitor<GastPti, GastPtiData>;

        GastPti();
        GastPti(const ModelDataT& data);
        ~GastPti() override;

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
                                GastPtiInternalVariables,
                                GastPtiExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        void    initModelParams(const ModelDataT& data);
        void    setDerivedParameters();
        void    initializeMonitor();
        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        RealT        R_{0.05};                    ///< Permanent droop
        RealT        T1_{0.4};                    ///< Fuel-valve time constant
        RealT        T2_{0.1};                    ///< Fuel-flow time constant
        RealT        T3_{3.0};                    ///< Exhaust-temperature time constant
        RealT        At_{1.0};                    ///< Ambient-temperature load limit
        RealT        Kt_{2.0};                    ///< Exhaust-temperature feedback gain
        RealT        Vmax_{1.0};                  ///< Maximum fuel-valve/turbine-power limit
        RealT        Vmin_{0.0};                  ///< Minimum fuel-valve/turbine-power limit
        RealT        Dturb_{0.0};                 ///< Turbine damping coefficient
        RealT        Trate_{100.0};               ///< Turbine-rating power base
        ResponseMode mode_{ResponseMode::Normal}; ///< Governor response mode

        // Derived parameters
        RealT va_component_base_{0};
        RealT sfix_{1.0}; ///< Fixed-mode dynamic selector

        ScalarT pref_set_{0};

        ComponentSignals<ScalarT, IdxT, GastPtiInternalVariables, GastPtiExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                           monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

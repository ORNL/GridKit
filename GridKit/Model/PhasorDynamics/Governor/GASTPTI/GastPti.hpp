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

      /// Indices into the GASTPTI state, derivative, and residual vectors.
      struct GastPtiIdx
      {
        static constexpr size_t XVALVE  = static_cast<size_t>(GastPtiInternalVariables::XVALVE);
        static constexpr size_t XFLOW   = static_cast<size_t>(GastPtiInternalVariables::XFLOW);
        static constexpr size_t XTEMP   = static_cast<size_t>(GastPtiInternalVariables::XTEMP);
        static constexpr size_t VLOAD   = static_cast<size_t>(GastPtiInternalVariables::VLOAD);
        static constexpr size_t VTEMP   = static_cast<size_t>(GastPtiInternalVariables::VTEMP);
        static constexpr size_t VLV     = static_cast<size_t>(GastPtiInternalVariables::VLV);
        static constexpr size_t PMECH   = static_cast<size_t>(GastPtiInternalVariables::PMECH);
        static constexpr size_t MAXIMUM = static_cast<size_t>(GastPtiInternalVariables::MAXIMUM);
      };

      /// Indices into the GASTPTI external-signal buffers.
      struct GastPtiExt
      {
        static constexpr size_t OMEGA   = static_cast<size_t>(GastPtiExternalVariables::OMEGA);
        static constexpr size_t PREF    = static_cast<size_t>(GastPtiExternalVariables::PREF);
        static constexpr size_t MAXIMUM = static_cast<size_t>(GastPtiExternalVariables::MAXIMUM);
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
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using ModelDataT = GastPtiData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<GastPti, GastPtiData>;

        GastPti();
        explicit GastPti(const ModelDataT& data);
        ~GastPti();

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
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        RealT        R_{static_cast<RealT>(0.05)};
        RealT        T1_{static_cast<RealT>(0.4)};
        RealT        T2_{static_cast<RealT>(0.1)};
        RealT        T3_{static_cast<RealT>(3.0)};
        RealT        At_{ONE<RealT>};
        RealT        Kt_{static_cast<RealT>(2.0)};
        RealT        Vmax_{ONE<RealT>};
        RealT        Vmin_{ZERO<RealT>};
        RealT        Dturb_{ZERO<RealT>};
        RealT        Trate_{static_cast<RealT>(100.0)};
        ResponseMode mode_{ResponseMode::Normal};

        RealT va_component_base_{ZERO<RealT>};
        RealT sfix_{ONE<RealT>};

        IdxT parameter_error_count_{0};

        ScalarT pref_set_{0};

        ComponentSignals<ScalarT, IdxT, GastPtiInternalVariables, GastPtiExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                           monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

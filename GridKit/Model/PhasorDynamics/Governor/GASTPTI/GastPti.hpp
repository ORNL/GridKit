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
#include <GridKit/Model/PhasorDynamics/Governor/ResponseMode.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /// Internal variables and residual rows of a `GastPti`.
      enum class GastPtiInternalVariables : size_t
      {
        XVALVE, ///< \f$x_V\f$ Differential fuel-valve state on component base [p.u.]
        XFLOW,  ///< \f$x_F\f$ Differential fuel-flow state on component base [p.u.]
        XTEMP,  ///< \f$x_T\f$ Differential exhaust-temperature feedback state on component base [p.u.]
        VLOAD,  ///< \f$V_D\f$ Algebraic speed/load fuel demand on component base [p.u.]
        VTEMP,  ///< \f$V_T\f$ Algebraic temperature-limit demand on component base [p.u.]
        VLV,    ///< \f$V\f$ Algebraic smooth low-value selector output on component base [p.u.]
        PMECH,  ///< \f$P_{\text{m}}\f$ Algebraic mechanical-power output on system base [p.u.]
        MAXIMUM ///< Number of GASTPTI internal variables and residual rows
      };

      /// External signal variables read or initialized by a `GastPti`.
      enum class GastPtiExternalVariables : size_t
      {
        OMEGA,  ///< \f$\omega\f$ Machine speed deviation [p.u.]
        PREF,   ///< \f$P^\mathrm{ref}\f$ Active-power/load reference on system base [p.u.]
        MAXIMUM ///< Number of GASTPTI external signal variables
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
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT         = GastPtiData<RealT, IdxT>;
        using MonitorT           = Model::VariableMonitor<GastPti, GastPtiData>;
        using InternalVariablesT = GastPtiInternalVariables;
        using ExternalVariablesT = GastPtiExternalVariables;

        GastPti();
        explicit GastPti(const ModelDataT& data);
        ~GastPti();

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
                                GastPtiInternalVariables,
                                GastPtiExternalVariables>&;

        const Model::VariableMonitorBase* getMonitor() const override;

        [[gnu::always_inline]] inline int evaluateInternalResidual(
            const ScalarT* y,
            const ScalarT* yp,
            const ScalarT* wb,
            const ScalarT* ws,
            ScalarT*       f);

      private:
        static constexpr size_t index(GastPtiInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        }

        static constexpr size_t index(GastPtiExternalVariables variable)
        {
          return static_cast<size_t>(variable);
        }

        static void checkConfiguration(bool condition, const char* message, int& errors);
        void        loadRealParameter(const ModelDataT& data,
                                      GastPtiParameters parameter,
                                      RealT&            target,
                                      const char*       name);
        bool        floorTimeConstant(RealT& value, const char* name);
        void        initializeParameters(const ModelDataT& data);
        void        initializeMonitor();
        void        setDerivedParameters();

        static RealT                              iramp(RealT value);
        [[gnu::always_inline]] inline scalar_type toComponentBase(scalar_type value) const;
        RealT                                     toSystemBase(RealT value) const;

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
        RealT Vmin_response_{ZERO<RealT>};
        RealT Vmax_response_{ONE<RealT>};
        RealT s_valve_{ONE<RealT>};

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

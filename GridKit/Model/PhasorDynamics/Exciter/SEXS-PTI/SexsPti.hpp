/**
 * @file SexsPti.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the SEXS-PTI exciter model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename real_type, typename index_type>
      struct SexsPtiData;
    } // namespace Exciter

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
    namespace Exciter
    {
      /// Internal variables of a `SexsPti`.
      enum class SexsPtiInternalVariables : size_t
      {
        VR,  ///< Lead-lag block state
        EFD, ///< Exciter field voltage output
        VTR, ///< Terminal voltage error signal
        MAXIMUM,
      };

      /// External variables of a `SexsPti`.
      enum class SexsPtiExternalVariables : size_t
      {
        VS, ///< Stabilizer output signal
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class SexsPti : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::time_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::J_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::residual_indices_;

      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using BusT       = BusBase<ScalarT, IdxT>;
        using ModelDataT = SexsPtiData<RealT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using MonitorT   = Model::VariableMonitor<SexsPti, SexsPtiData>;

        SexsPti(BusT* bus);
        SexsPti(BusT* bus, const ModelDataT& data);
        ~SexsPti();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                SexsPtiInternalVariables,
                                SexsPtiExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        FORCE_INLINE int evaluateInternalResidual(
            ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        BusT* bus_{nullptr};

        RealT Ta_{0};
        RealT Tb_{0};
        RealT Te_{0};
        RealT K_{0};
        RealT Efdmax_{0};
        RealT Efdmin_{0};

        int missing_param_count_{0};

        ScalarT vref_{0};
        ScalarT vOEL_{0};
        ScalarT vUEL_{0};

        ComponentSignals<ScalarT, IdxT, SexsPtiInternalVariables, SexsPtiExternalVariables> signals_;

        std::unique_ptr<MonitorT> monitor_;

        void initModelParams(const ModelDataT& data);
        void initializeMonitor();

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

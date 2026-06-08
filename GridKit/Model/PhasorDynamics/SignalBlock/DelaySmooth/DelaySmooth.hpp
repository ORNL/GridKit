/**
 * @file DelaySmooth.hpp
 * @brief Declaration of the smooth lag-chain delay signal block.
 */

#pragma once

#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class SignalNode;

    namespace SignalBlock
    {
      template <typename real_type, typename index_type>
      struct DelaySmoothData;

      enum class DelaySmoothInternalVariables : size_t
      {
        OUT,
        MAXIMUM,
      };

      enum class DelaySmoothExternalVariables : size_t
      {
        IN,
        MAXIMUM,
      };

      /// Smooth delay signal block: approximates a transport delay tau by a
      /// string of N first-order lag blocks. The block owns N differential
      /// states, has a hand-written banded Jacobian, and keeps no history.
      template <typename scalar_type, typename index_type>
      class DelaySmooth : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::J_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::time_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT = DelaySmoothData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<DelaySmooth, DelaySmoothData>;

        DelaySmooth();
        DelaySmooth(const ModelDataT& data);
        ~DelaySmooth();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT, IdxT, DelaySmoothInternalVariables, DelaySmoothExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        int evaluateInternalResidual(ScalarT* y, ScalarT* yp, ScalarT* ws, ScalarT* f);

      private:
        void  initializeParameters(const ModelDataT& data);
        void  initializeMonitor();
        RealT inputValue() const;

        RealT delay_{0.0};
        RealT dt_min_{0.0};
        RealT T_{0.0};
        RealT G_{0.0};
        IdxT  N_{0};

        std::vector<ScalarT> ws_;

        ComponentSignals<ScalarT, IdxT, DelaySmoothInternalVariables, DelaySmoothExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                                   monitor_;
      };
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit

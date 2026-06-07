/**
 * @file Delay.hpp
 * @brief Declaration of the transport-delay signal block.
 */

#pragma once

#include <memory>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/SignalBlock/Delay/SignalHistory.hpp>
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
      struct DelayData;

      enum class DelayInternalVariables : size_t
      {
        OUT,
        MAXIMUM,
      };

      enum class DelayExternalVariables : size_t
      {
        IN,
        MAXIMUM,
      };

      /// Transport-delay signal block: writes its input signal, delayed by a
      /// fixed lag tau, to its output signal. The block owns no solver state and
      /// has no residual equations; it records its input at accepted steps and
      /// emits the interpolated past value.
      template <typename scalar_type, typename index_type>
      class Delay : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::nnz_;
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
        using ModelDataT = DelayData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Delay, DelayData>;

        Delay();
        Delay(const ModelDataT& data);
        ~Delay();

        int  setGridKitComponentID(IdxT component_id) override final;
        int  allocate() override final;
        int  verify() const override final;
        int  initialize() override final;
        int  tagDifferentiable() override final;
        int  evaluateResidual() override final;
        int  evaluateJacobian() override final;
        void updateTime(RealT t, RealT a) override final;

        // Step-history lifecycle hooks (see Evaluator/Component).
        int   resetHistory(RealT t0) override final;
        int   stepAccepted(RealT t) override final;
        RealT maxStepSize() const override final;

        auto getSignals()
            -> ComponentSignals<ScalarT, IdxT, DelayInternalVariables, DelayExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

      private:
        void  initializeParameters(const ModelDataT& data);
        void  initializeMonitor();
        RealT inputValue() const;

        RealT   delay_{0.0};
        RealT   prehistory_value_{0.0};
        bool    has_prehistory_{false};
        ScalarT value_{0.0};
        IdxT    invalid_index_{INVALID_INDEX<IdxT>};

        SignalHistory<RealT> history_;

        ComponentSignals<ScalarT, IdxT, DelayInternalVariables, DelayExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;
      };
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit

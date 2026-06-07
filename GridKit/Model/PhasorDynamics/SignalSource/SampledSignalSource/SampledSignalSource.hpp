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

    namespace SignalSource
    {
      template <typename real_type, typename index_type>
      struct SampledSignalSourceData;

      enum class SampledSignalSourceInternalVariables : size_t
      {
        OUT,
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class SampledSignalSource : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::J_;
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
        using ModelDataT = SampledSignalSourceData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<SampledSignalSource, SampledSignalSourceData>;

        SampledSignalSource();
        SampledSignalSource(const ModelDataT& data);
        ~SampledSignalSource();

        int  setGridKitComponentID(IdxT component_id) override final;
        int  allocate() override final;
        int  verify() const override final;
        int  initialize() override final;
        int  tagDifferentiable() override final;
        int  evaluateResidual() override final;
        int  evaluateJacobian() override final;
        void updateTime(RealT t, RealT a) override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                SampledSignalSourceInternalVariables,
                                NoVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        RealT sample(RealT t) const;
        RealT valueAt(RealT t) const;

      private:
        void initializeParameters(const ModelDataT& data);
        void loadSamples(const ModelDataT& data);
        void loadCsvSamples(const ModelDataT& data);
        void validateSamples();
        void initializeMonitor();

        RealT   scale_{1.0};
        RealT   offset_{0.0};
        ScalarT value_{0.0};
        IdxT    invalid_index_{INVALID_INDEX<IdxT>};
        int     configuration_errors_{0};

        std::vector<std::pair<RealT, RealT>> samples_;

        ComponentSignals<ScalarT, IdxT, SampledSignalSourceInternalVariables, NoVariables> signals_;
        std::unique_ptr<MonitorT>                                                          monitor_;
      };
    } // namespace SignalSource
  } // namespace PhasorDynamics
} // namespace GridKit

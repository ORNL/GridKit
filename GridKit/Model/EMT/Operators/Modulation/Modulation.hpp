#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/ModulationData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// DC normalization of three-phase voltage commands without DAE variables.
    template <typename scalar_type, typename index_type>
    class Modulation : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using SignalT    = Signal<ScalarT, IdxT>;
      using ModelDataT = ModulationData<RealT, IdxT>;
      using Outputs    = typename ModelDataT::Outputs;
      using MonitorT   = Model::VariableMonitor<Modulation, ModulationData>;

      Modulation();
      explicit Modulation(const ModelDataT& data);
      ~Modulation() override;

      int setGridKitComponentID(IdxT id) override final;
      int allocate() override final;
      int verify() const override final;

      int initialize(const std::map<Outputs, RealT>& outputs = {});

      typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override
      {
        return {};
      }

      int initializeState(const std::map<std::string, RealT>& values) override
      {
        return this->initializeOutputs(*this, values);
      }

      void validateInitialState(const std::map<std::string, RealT>& values) const override
      {
        this->template parseInitialOutputs<Modulation>(values);
      }

      int setAbsoluteTolerance(RealT) override final;
      int evaluateInternalResidual() override final;
      int evaluateExternalResidual() override final;
      int evaluateResidual() override final;
      int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      /// Publish a modulation output on a named scalar signal. No DAE index is assigned.
      void    assignOutput(Outputs output, SignalT* signal);
      ScalarT output(Outputs output) const;
      void    attachInput(const std::array<SignalT*, 3>& command, SignalT* vdc);

      SignalT& outputSignal(Outputs output)
      {
        return output_port_.at(static_cast<size_t>(output));
      }

    private:
      void                              appendOutputGradient(Outputs output, typename SignalT::GradientT& gradient, RealT scale) const;
      const Model::VariableMonitorBase* getMonitor() const override;

      std::array<SignalT*, 4>   input_{};
      std::array<SignalT, 3>    output_port_;
      std::array<SignalT*, 3>   assigned_output_{};
      std::unique_ptr<MonitorT> monitor_;
    };
  } // namespace EMT
} // namespace GridKit

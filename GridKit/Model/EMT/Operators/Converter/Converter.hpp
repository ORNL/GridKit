#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Converter/ConverterData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Two-level voltage-source bridge with no internal variables or residual rows.
    template <typename scalar_type, typename index_type>
    class Converter : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using SignalT    = Signal<ScalarT, IdxT>;
      using ModelDataT = ConverterData<RealT, IdxT>;
      using Outputs    = typename ModelDataT::Outputs;
      using MonitorT   = Model::VariableMonitor<Converter, ConverterData>;

      Converter();
      explicit Converter(const ModelDataT& data);
      ~Converter() override;

      int setGridKitComponentID(IdxT id) override final;
      int allocate() override final;
      int verify() const override final;

      int initialize(const std::map<Outputs, RealT>& outputs = {});

      int initializeState(const std::map<std::string, RealT>& values) override
      {
        return this->initializeOutputs(*this, values);
      }

      int initializationOrder() const noexcept override final
      {
        return 4;
      }

      int setAbsoluteTolerance(RealT) override final;
      int evaluateInternalResidual() override final;
      int evaluateExternalResidual() override final;
      int evaluateResidual() override final;
      int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      /// Publish one phase on a named scalar signal. No DAE index is assigned.
      void    assignOutput(size_t phase, SignalT* signal);
      ScalarT output(size_t phase) const;
      void    attachInput(SignalT* a, SignalT* b, SignalT* c, SignalT* vdc);

      SignalT& outputSignal(Outputs output)
      {
        return output_port_.at(static_cast<size_t>(output));
      }

      /// Two-level bridge projection, also usable directly in residual kernels.
      __attribute__((always_inline)) static ABCVector<ScalarT> voltage(const ABCVector<ScalarT>& s, ScalarT vdc)
      {
        return {vdc * ((s[0] - s[1]) + (s[0] - s[2])) / 3,
                vdc * ((s[1] - s[0]) + (s[1] - s[2])) / 3,
                vdc * ((s[2] - s[0]) + (s[2] - s[1])) / 3};
      }

    private:
      void                              appendOutputGradient(size_t phase, typename SignalT::GradientT& gradient, RealT scale) const;
      const Model::VariableMonitorBase* getMonitor() const override;

      std::array<SignalT*, 4>   input_{};
      std::array<SignalT, 3>    output_port_;
      std::array<SignalT*, 3>   assigned_output_{};
      std::unique_ptr<MonitorT> monitor_;
    };
  } // namespace EMT
} // namespace GridKit

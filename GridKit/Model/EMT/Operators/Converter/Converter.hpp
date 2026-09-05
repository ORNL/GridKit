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
      using Port3T     = Port3<ScalarT, IdxT>;
      using ModelDataT = ConverterData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<Converter, ConverterData>;

      Converter();
      explicit Converter(const ModelDataT& data);
      ~Converter() override;

      int setGridKitComponentID(IdxT id) override final;
      int allocate() override final;
      int verify() const override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int setAbsoluteTolerance(RealT) override final;
      int evaluateInternalResidual() override final;
      int evaluateExternalResidual() override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      /// Publish one phase on a named scalar signal. No DAE index is assigned.
      void    assignOutput(size_t phase, SignalT* signal);
      ScalarT output(size_t phase) const;
      void    attachInput(SignalT* a, SignalT* b, SignalT* c, SignalT* vdc);
      void    attachInput(Port3T* s, SignalT* vdc);

      Port3T& voltagePort()
      {
        return output_port_;
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
      Port3T                    output_port_;
      std::array<SignalT*, 3>   assigned_output_{};
      std::unique_ptr<MonitorT> monitor_;
    };
  } // namespace EMT
} // namespace GridKit

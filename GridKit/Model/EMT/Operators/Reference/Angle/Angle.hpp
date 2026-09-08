#pragma once

#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Angle/AngleData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Electrical reference angle with one owned differential state.
    template <typename scalar_type, typename index_type>
    class Angle : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using Base       = Component<ScalarT, IdxT>;
      using RealT      = typename Base::RealT;
      using SignalT    = typename Base::SignalT;
      using ModelDataT = AngleData<RealT, IdxT>;
      using Outputs    = typename ModelDataT::Outputs;
      using MonitorT   = Model::VariableMonitor<Angle, AngleData>;

      explicit Angle(const ModelDataT& data);
      ~Angle() override;

      int setGridKitComponentID(IdxT id) override final;
      int allocate() override final;
      int verify() const override final;

      int                                 initialize(const std::map<Outputs, RealT>& outputs = {});
      void                                validateInitialState(const std::map<std::string, RealT>& values) const override;
      int                                 initializeState(const std::map<std::string, RealT>& values) override;
      typename Base::InitializationPortsT initializationPorts() override;
      void                                prepareInitialization(typename Base::InitialStateT& initial) override;

      int setAbsoluteTolerance(RealT tolerance) override final;
      int evaluateInternalResidual() override final;
      int evaluateResidual() override final;
      int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      void     attachInput(SignalT* omega);
      void     assignOutput(Outputs output, SignalT* signal);
      SignalT& outputSignal(Outputs output);

    private:
      const Model::VariableMonitorBase* getMonitor() const override;

      SignalT                   angle_;
      SignalT*                  alias_{nullptr};
      SignalT*                  omega_{nullptr};
      size_t                    capacity_{0};
      std::unique_ptr<MonitorT> monitor_;
    };
  } // namespace EMT
} // namespace GridKit

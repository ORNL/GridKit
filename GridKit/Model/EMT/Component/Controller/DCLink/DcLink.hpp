#pragma once

#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/DCLink/DcLinkData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Ideal DC-link capacitor with one owned voltage state.
      template <typename scalar_type, typename index_type>
      class DcLink : public Component<scalar_type, index_type>
      {
      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using Base       = Component<ScalarT, IdxT>;
        using RealT      = typename Base::RealT;
        using SignalT    = typename Base::SignalT;
        using ModelDataT = DcLinkData<RealT, IdxT>;
        using Outputs    = typename ModelDataT::Outputs;
        using MonitorT   = Model::VariableMonitor<DcLink, DcLinkData>;

        explicit DcLink(const ModelDataT& data);
        ~DcLink() override;

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

        /// Currents are positive into the link from the source and out to the converter.
        void     attachInput(SignalT* isrc, SignalT* idc);
        void     assignOutput(Outputs output, SignalT* signal);
        SignalT& outputSignal(Outputs output);

      private:
        const Model::VariableMonitorBase* getMonitor() const override;

        RealT                     capacitance_;
        SignalT                   voltage_;
        SignalT*                  alias_{nullptr};
        SignalT*                  isrc_{nullptr};
        SignalT*                  idc_{nullptr};
        size_t                    capacity_{0};
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

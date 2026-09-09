/**
 * @file Bus.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT bus model.
 */
#pragma once

#include <GridKit/Model/EMT/Component/Bus/KCL.hpp>
#include <GridKit/Model/EMT/Component/Source/Norton/Norton.hpp>
#include <GridKit/Model/EMT/Container.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class Bus : public Container<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using Base         = Container<ScalarT, IdxT>;
      using RealT        = typename Base::RealT;
      using SignalT      = typename Base::SignalT;
      using ModelDataT   = BusData<RealT, IdxT>;
      using Outputs      = typename ModelDataT::Outputs;
      using KCLT         = KCL<ScalarT, IdxT>;
      using NortonT      = Norton<ScalarT, IdxT>;
      using AdmittanceT  = VectorFit<ScalarT, IdxT>;
      using PhaseSignals = typename NortonT::PhaseSignals;
      using PhaseOrder   = typename KCLT::PhaseOrder;
      using YDataT       = typename NortonT::YDataT;
      using MonitorT     = Model::VariableMonitor<Bus, BusData>;

      Bus();
      explicit Bus(const ModelDataT& data);
      ~Bus() override;

      int initialize(const std::map<Outputs, RealT>& outputs = {});

      int initialize(const std::map<std::string, std::map<std::string, RealT>>& state, RealT omega = RealT{0}) override
      {
        return Base::initialize(state, omega);
      }

      int initializeSteadyState(RealT omega);
      int setAbsoluteTolerance(RealT tolerance) override;

      void addCurrent(size_t phase, SignalT& signal, RealT sign = ONE<RealT>)
      {
        kcl_.addCurrent(phase, signal, sign);
      }

      void assignOutput(Outputs output, SignalT* signal)
      {
        kcl_.assignOutput(output, signal);
      }

      using Base::outputSignal;

      SignalT& outputSignal(Outputs output)
      {
        return kcl_.outputSignal(output);
      }

      PhaseSignals voltages(PhaseOrder phases = {0, 1, 2})
      {
        return kcl_.voltages(phases);
      }

      IdxT voltagePhase(const SignalT* signal) const
      {
        return kcl_.voltagePhase(signal);
      }

      NortonT&     addNorton(std::string name, const YDataT& Y, PhaseSignals incident = {}, RealT scale = ONE<RealT>, PhaseOrder phases = {0, 1, 2});
      AdmittanceT& addShunt(std::string name, const YDataT& Y, RealT scale = ONE<RealT>, PhaseOrder phases = {0, 1, 2});

      NortonT& norton(std::string_view name)
      {
        return this->template component<NortonT>(name);
      }

    protected:
      typename Base::ComponentT* initialStateComponent() override
      {
        return &kcl_;
      }

    private:
      void                              initializeMonitor();
      const Model::VariableMonitorBase* getMonitor() const override;

      RealT                                                     v_scale_{ONE<RealT>};
      KCLT&                                                     kcl_;
      SignalT                                                   zero_;
      std::array<std::vector<SignalT*>, 3>                      shunt_monitors_;
      std::array<std::vector<std::pair<AdmittanceT*, IdxT>>, 3> shunts_;
      std::unique_ptr<MonitorT>                                 monitor_;
    };
  } // namespace EMT
} // namespace GridKit

/**
 * @file Filter.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Three-phase LCL filter between converter voltage signals and an EMT bus.
 */
#pragma once

#include <array>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Filter/FilterData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class FilterInternalVariables : size_t
    {
      IA,
      IB,
      IC,
      VOA,
      VOB,
      VOC,
      IGA,
      IGB,
      IGC,
      MAXIMUM,
    };

    enum class FilterExternalVariables : size_t
    {
      VA,
      VB,
      VC,
      EA,
      EB,
      EC,
      MAXIMUM,
    };

    template <typename scalar_type, typename index_type>
    class Filter : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using Base         = Component<ScalarT, IdxT>;
      using RealT        = typename Base::RealT;
      using ModelDataT   = FilterData<RealT, IdxT>;
      using Outputs      = typename ModelDataT::Outputs;
      using SignalT      = typename Base::SignalT;
      using MonitorT     = Model::VariableMonitor<Filter, FilterData>;
      using PhaseSignals = std::array<SignalT*, 3>;

      explicit Filter(const ModelDataT& data);
      ~Filter() override;

      void     attachInput(PhaseSignals voltage, PhaseSignals source);
      void     assignOutput(Outputs output, SignalT* signal);
      SignalT& outputSignal(Outputs output);

      SignalT& inputSignal(FilterInputs input)
      {
        return *signals_.getAttachedSignal(static_cast<FilterExternalVariables>(input));
      }

      /// Grid-side current, positive into the connected bus.
      SignalT& currentSignal(size_t phase)
      {
        return output_.at(static_cast<size_t>(Outputs::iga) + phase);
      }

      int                                 setGridKitComponentID(IdxT id) override final;
      int                                 allocate() override final;
      int                                 verify() const override final;
      int                                 initialize(const std::map<Outputs, RealT>& outputs = {});
      int                                 initializeState(const std::map<std::string, RealT>& values) override;
      void                                validateInitialState(const std::map<std::string, RealT>& values) const override;
      typename Base::InitializationPortsT initializationPorts() override;
      void                                prepareInitialization(typename Base::InitialStateT& initial) override;
      int                                 setAbsoluteTolerance(RealT tolerance) override final;
      int                                 evaluateInternalResidual() override final;
      int                                 evaluateResidual() override final;
      int                                 assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      auto& getSignals()
      {
        return signals_;
      }

      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      void                              initializeParameters(const ModelDataT& data);
      void                              initializeMonitor();
      const Model::VariableMonitorBase* getMonitor() const override;

      ABCMatrix<RealT> Rs_{};
      ABCMatrix<RealT> Ls_{};
      ABCMatrix<RealT> C_{};
      ABCMatrix<RealT> Rg_{};
      ABCMatrix<RealT> Lg_{};

      ComponentSignals<ScalarT, IdxT, FilterInternalVariables, FilterExternalVariables> signals_;
      std::array<SignalT, 9>                                                            output_;
      std::array<SignalT*, 9>                                                           alias_{};
      std::unique_ptr<MonitorT>                                                         monitor_;
    };
  } // namespace EMT
} // namespace GridKit

/**
 * @file Pll.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Synchronous-reference-frame phase-locked loop.
 */
#pragma once

#include <array>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/EMT/Operators/Reference/PLL/PllData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class PllInternalVariables : size_t
    {
      THETA,
      XI,
      OMEGA,
      MAXIMUM
    };
    enum class PllExternalVariables : size_t
    {
      VA,
      VB,
      VC,
      MAXIMUM
    };

    template <typename scalar_type, typename index_type>
    class Pll : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::time_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::y_ext_;
      using Component<scalar_type, index_type>::yp_ext_;
      using Component<scalar_type, index_type>::variable_indices_ext_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::allocated_;

    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT   = PllData<RealT, IdxT>;
      using Outputs      = typename ModelDataT::Outputs;
      using SignalT      = Signal<ScalarT, IdxT>;
      using MonitorT     = Model::VariableMonitor<Pll, PllData>;
      using PhaseSignals = std::array<SignalT*, 3>;

      Pll();
      explicit Pll(const ModelDataT& data);
      ~Pll();
      void     attachInput(PhaseSignals voltage);
      void     assignOutput(Outputs output, SignalT* signal);
      SignalT& outputSignal(Outputs output);

      SignalT& inputSignal(PllInputs input)
      {
        return *signals_.getAttachedSignal(static_cast<PllExternalVariables>(input));
      }

      int                                                     setGridKitComponentID(IdxT) override final;
      int                                                     allocate() override final;
      int                                                     verify() const override final;
      int                                                     initialize(const std::map<Outputs, RealT>& outputs = {});
      int                                                     initializeState(const std::map<std::string, RealT>& values) override;
      void                                                    validateInitialState(const std::map<std::string, RealT>& values) const override;
      typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override;
      int                                                     setAbsoluteTolerance(RealT) override final;
      int                                                     evaluateInternalResidual() override final;
      int                                                     evaluateResidual() override final;
      int                                                     assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      auto& getSignals()
      {
        return signals_;
      }

      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      void                                          initializeParameters(const ModelDataT& data);
      void                                          initializeMonitor();
      const Model::VariableMonitorBase*             getMonitor() const override;
      __attribute__((always_inline)) inline ScalarT quadrature(ScalarT theta, const ScalarT* voltage) const;

      RealT                                                                       V_{0.0};
      RealT                                                                       freq_{0.0};
      RealT                                                                       Kp_{0.0};
      RealT                                                                       Ki_{0.0};
      RealT                                                                       omega0_{0.0};
      RealT                                                                       projection_scale_{0.0};
      RealT                                                                       beta_scale_{0.0};
      ComponentSignals<ScalarT, IdxT, PllInternalVariables, PllExternalVariables> signals_;
      std::array<SignalT, 2>                                                      output_;
      std::array<SignalT*, 2>                                                     alias_{};
      std::unique_ptr<MonitorT>                                                   monitor_;
    };
  } // namespace EMT
} // namespace GridKit

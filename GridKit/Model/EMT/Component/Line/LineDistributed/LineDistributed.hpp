/**
 * @file LineDistributed.hpp
 * @brief Declaration of the EMT distributed line model.
 */
#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Line/LineDistributed/LineDistributedData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Propagation/Propagation.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `LineDistributed`.
    enum class LineDistributedInternalVariables : size_t
    {
      IREF1A, ///< Reflected current at terminal 1, phase a
      IREF1B, ///< Reflected current at terminal 1, phase b
      IREF1C, ///< Reflected current at terminal 1, phase c
      IREF2A, ///< Reflected current at terminal 2, phase a
      IREF2B, ///< Reflected current at terminal 2, phase b
      IREF2C, ///< Reflected current at terminal 2, phase c
      MAXIMUM,
    };

    /// External variables of a `LineDistributed` equation block.
    enum class LineDistributedExternalVariables : size_t
    {
      IC1A,   ///< Characteristic-admittance current at terminal 1, phase a
      IC1B,   ///< Characteristic-admittance current at terminal 1, phase b
      IC1C,   ///< Characteristic-admittance current at terminal 1, phase c
      IC2A,   ///< Characteristic-admittance current at terminal 2, phase a
      IC2B,   ///< Characteristic-admittance current at terminal 2, phase b
      IC2C,   ///< Characteristic-admittance current at terminal 2, phase c
      IINC1A, ///< Propagation output at terminal 1, phase a
      IINC1B, ///< Propagation output at terminal 1, phase b
      IINC1C, ///< Propagation output at terminal 1, phase c
      IINC2A, ///< Propagation output at terminal 2, phase a
      IINC2B, ///< Propagation output at terminal 2, phase b
      IINC2C, ///< Propagation output at terminal 2, phase c
      MAXIMUM,
    };

    /**
     * @brief Reciprocal three-phase distributed EMT line.
     *
     * The buses own the characteristic-admittance currents and states.
     * The line owns the reflected currents and independent propagation
     * states and accepted histories in each direction.
     */
    template <typename scalar_type, typename index_type>
    class LineDistributed : public Component<scalar_type, index_type>
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
      using Component<scalar_type, index_type>::equation_size_;

    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT   = LineDistributedData<RealT, IdxT>;
      using Outputs      = typename ModelDataT::Outputs;
      using SignalT      = Signal<ScalarT, IdxT>;
      using PropagationT = Propagation<ScalarT, IdxT>;
      using MonitorT     = Model::VariableMonitor<LineDistributed, LineDistributedData>;
      using PhaseSignals = std::array<SignalT*, 3>;

      LineDistributed();
      explicit LineDistributed(const ModelDataT& data);
      ~LineDistributed() override;

      void     attachTerminal(size_t end, PhaseSignals characteristic);
      SignalT& incidentSignal(size_t end, size_t phase);
      SignalT& outputSignal(Outputs output);
      void     assignOutput(Outputs output, SignalT* signal);

      int  setGridKitComponentID(IdxT component_id) override final;
      int  allocate() override final;
      int  verify() const override final;
      int  initialize(const std::map<Outputs, RealT>& outputs = {});
      void setPrehistory(RealT omega, const std::array<ABCVector<RealT>, 2>& value, const std::array<ABCVector<RealT>, 2>& derivative);

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
        const auto outputs = this->template parseInitialOutputs<LineDistributed>(values);
        for (const auto& [output, value] : outputs)
          if (static_cast<size_t>(output) >= static_cast<size_t>(LineDistributedInternalVariables::MAXIMUM))
            throw std::invalid_argument("LineDistributed initial outputs must be reflected currents");
      }

      int setAbsoluteTolerance(RealT rel_tol) override final;
      int evaluateInternalResidual() override final;
      int evaluateResidual() override final;
      int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      auto getSignals() -> ComponentSignals<ScalarT,
                                            IdxT,
                                            LineDistributedInternalVariables,
                                            LineDistributedExternalVariables>&
      {
        return signals_;
      }

    private:
      void                              initializeParameters(const ModelDataT& data);
      void                              initializePorts();
      void                              initializeMonitor();
      void                              setDerivedParams();
      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Input parameters */
      IdxT                            N_{3};
      IdxT                            K_{3};
      ABCVector<IdxT>                 conductors_{{1, 2, 3}};
      RealT                           omega_{0.0};
      std::array<ABCVector<RealT>, 2> history_{};
      std::array<ABCVector<RealT>, 2> history_derivative_{};
      bool                            has_history_{false};

      /* Propagation operators and reflected-current ports */
      std::array<std::unique_ptr<PropagationT>, 2>                                                        propagation_;
      std::array<std::array<SignalT, 3>, 2>                                                               reflected_;
      ComponentSignals<ScalarT, IdxT, LineDistributedInternalVariables, LineDistributedExternalVariables> signals_;
      std::unique_ptr<MonitorT>                                                                           monitor_;
      size_t                                                                                              jacobian_capacity_{0};
    };

  } // namespace EMT
} // namespace GridKit

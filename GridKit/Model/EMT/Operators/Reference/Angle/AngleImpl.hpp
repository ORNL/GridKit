#pragma once

#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Angle/Angle.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Angle<scalar_type, index_type>::Angle(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      this->equation_size_ = this->size_ = 1;
      angle_.claimProducer();
      using Mon = typename ModelDataT::MonitorableVariables;
      monitor_->set(Mon::theta, [this]
                    { return angle_.read(); });
    }

    template <typename scalar_type, typename index_type>
    typename Angle<scalar_type, index_type>::SignalT& Angle<scalar_type, index_type>::outputSignal(Outputs output)
    {
      if (output != Outputs::theta)
        throw std::invalid_argument("Angle: invalid output");
      return angle_;
    }

    template <typename scalar_type, typename index_type>
    void Angle<scalar_type, index_type>::attachInput(SignalT* omega)
    {
      if (this->allocated_)
        throw std::logic_error("Angle inputs cannot change after allocation");
      if (!omega)
        throw std::invalid_argument("Angle requires omega");
      omega_ = omega;
    }

    template <typename scalar_type, typename index_type>
    void Angle<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      if (this->allocated_)
        throw std::logic_error("Angle outputs cannot change after allocation");
      if (output != Outputs::theta || !signal || alias_)
        throw std::invalid_argument("Angle: invalid or duplicate output assignment");
      signal->claimProducer();
      alias_ = signal;
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      this->gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::allocate()
    {
      if (!omega_)
        throw std::invalid_argument("Angle requires omega");
      if (!this->allocated_)
        this->allocateVectors(1);
      this->tag_.resize(1);
      this->variable_indices_.resize(1);
      this->residual_indices_.resize(1);
      this->assignGlobalIndices(0);
      this->bindSignal(angle_, 0);
      if (alias_)
        this->bindSignal(*alias_, 0);
      this->allocateExternalVectors(1, 0);
      this->setExternalVariableSignal(0, omega_);
      this->allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::verify() const
    {
      return !omega_ || !omega_->linked();
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      angle_.init(static_cast<ScalarT>(this->outputValue(outputs, Outputs::theta, RealT{0})));
      angle_.initDerivative(ScalarT{0});
      this->y_.setDataUpdated();
      this->yp_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Angle<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
    {
      this->template parseInitialOutputs<Angle>(values);
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
    {
      return this->initializeOutputs(*this, values);
    }

    template <typename scalar_type, typename index_type>
    typename Component<scalar_type, index_type>::InitializationPortsT Angle<scalar_type, index_type>::initializationPorts()
    {
      return {{}, {{"theta", &angle_}}, {}};
    }

    template <typename scalar_type, typename index_type>
    void Angle<scalar_type, index_type>::prepareInitialization(typename Base::InitialStateT& initial)
    {
      const auto outputs = this->template parseInitialOutputs<Angle>(initial.outputs(*this));
      const auto value   = this->outputValue(outputs, Outputs::theta, RealT{0});
      initial.provide(angle_, value);
      if (alias_)
        initial.provide(*alias_, value);
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->f_.getData()[0] = -this->yp_.getData()[0] + omega_->read();
      this->f_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    int Angle<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      typename SignalT::GradientT gradient;
      if (y_scale != ZERO<RealT>)
        omega_->appendGradient(gradient, y_scale);
      if (yp_scale != ZERO<RealT>)
        gradient.emplace_back(this->getVariableIndex(0), -yp_scale);
      const auto entries = gradient.size();
      if (entries != capacity_)
      {
        this->resetJacobianStructure();
        delete[] this->J_rows_buffer_;
        delete[] this->J_cols_buffer_;
        delete[] this->J_vals_buffer_;
        this->J_rows_buffer_ = new IdxT[entries];
        this->J_cols_buffer_ = new IdxT[entries];
        this->J_vals_buffer_ = new RealT[entries];
        capacity_            = entries;
      }
      this->nnz_ = 0;
      for (const auto& [column, value] : gradient)
      {
        const auto j            = this->nnz_++;
        this->J_rows_buffer_[j] = this->getResidualIndex(0);
        this->J_cols_buffer_[j] = column;
        this->J_vals_buffer_[j] = value;
      }
      return entries == 0 ? 0 : this->constructCoo();
    }

    template <typename scalar_type, typename index_type>
    Angle<scalar_type, index_type>::~Angle() = default;

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Angle<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }
  } // namespace EMT
} // namespace GridKit

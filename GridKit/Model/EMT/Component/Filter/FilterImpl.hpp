#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include <GridKit/Model/EMT/Component/Filter/Filter.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Filter<scalar_type, index_type>::Filter(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      this->equation_size_ = this->size_ = 9;
      for (size_t n = 0; n < output_.size(); ++n)
        signals_.assignSignal(static_cast<FilterInternalVariables>(n), &output_[n]);
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Filter<scalar_type, index_type>::~Filter() = default;

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using P = typename ModelDataT::Parameters;
      Rs_     = parameter<ABCMatrix<RealT>>(data, P::Rs, Rs_);
      Ls_     = parameter<ABCMatrix<RealT>>(data, P::Ls);
      C_      = parameter<ABCMatrix<RealT>>(data, P::C);
      Rg_     = parameter<ABCMatrix<RealT>>(data, P::Rg, Rg_);
      Lg_     = parameter<ABCMatrix<RealT>>(data, P::Lg);

      // Principal minors characterize positive (semi)definite symmetric
      // three-phase matrices. Normalize to keep the check independent of units.
      auto passive = [](const ABCMatrix<RealT>& matrix, bool strict)
      {
        RealT scale = 0;
        for (const auto& row : matrix)
          for (const auto value : row)
          {
            if (!std::isfinite(value))
              return false;
            scale = std::max(scale, std::abs(value));
          }
        if (scale == RealT{0})
          return !strict;
        ABCMatrix<RealT> a{};
        const RealT      tolerance = 3 * std::numeric_limits<RealT>::epsilon();
        for (size_t n = 0; n < 3; ++n)
          for (size_t k = 0; k < 3; ++k)
          {
            a[n][k] = matrix[n][k] / scale;
            if (std::abs(a[n][k] - matrix[k][n] / scale) > tolerance)
              return false;
          }
        auto positive = [=](RealT value)
        { return strict ? value > tolerance : value >= -tolerance; };
        for (size_t n = 0; n < 3; ++n)
        {
          if (!positive(a[n][n]))
            return false;
          for (size_t k = n + 1; k < 3; ++k)
            if (!positive(a[n][n] * a[k][k] - a[n][k] * a[k][n]))
              return false;
        }
        return positive(a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
      };
      if (!passive(Rs_, false) || !passive(Rg_, false))
        throw std::invalid_argument("Filter: Rs and Rg must be finite symmetric positive semidefinite matrices");
      if (!passive(Ls_, true) || !passive(C_, true) || !passive(Lg_, true))
        throw std::invalid_argument("Filter: Ls, C and Lg must be finite symmetric positive definite matrices");
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::attachInput(PhaseSignals voltage, PhaseSignals source)
    {
      if (this->allocated_)
        throw std::logic_error("Filter inputs cannot change after allocation");
      for (size_t p = 0; p < 3; ++p)
      {
        signals_.attachSignal(static_cast<FilterExternalVariables>(p), voltage[p]);
        signals_.attachSignal(static_cast<FilterExternalVariables>(3 + p), source[p]);
      }
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      if (this->allocated_)
        throw std::logic_error("Filter outputs cannot change after allocation");
      const auto n = static_cast<size_t>(output);
      if (n >= output_.size() || !signal || alias_[n])
        throw std::invalid_argument("Filter: invalid or duplicate output assignment");
      signal->claimProducer();
      alias_[n] = signal;
    }

    template <typename scalar_type, typename index_type>
    typename Filter<scalar_type, index_type>::SignalT& Filter<scalar_type, index_type>::outputSignal(Outputs output)
    {
      return output_.at(static_cast<size_t>(output));
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      this->gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::allocate()
    {
      using V = FilterExternalVariables;
      if (signals_.attachedSignals({V::VA, V::VB, V::VC, V::EA, V::EB, V::EC}).size() != 6)
        throw std::invalid_argument("Filter: all voltage inputs are required");
      if (!this->allocated_)
        this->allocateVectors(this->size_);
      this->tag_.resize(9);
      this->variable_indices_.resize(9);
      this->residual_indices_.resize(9);
      this->assignGlobalIndices(0);
      this->allocateExternalVectors(6, 0);
      signals_.registerExternalVariableSignals(*this);
      signals_.bindInternalVariableSignals(*this);
      for (size_t n = 0; n < alias_.size(); ++n)
        if (alias_[n])
          this->bindSignal(*alias_[n], static_cast<IdxT>(n));
      this->allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::verify() const
    {
      using V           = FilterExternalVariables;
      const auto inputs = signals_.attachedSignals({V::VA, V::VB, V::VC, V::EA, V::EB, V::EC});
      int        errors = inputs.size() != 6;
      for (const auto* input : inputs)
        errors += !input->linked();
      return errors;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      for (size_t n = 0; n < output_.size(); ++n)
      {
        output_[n].init(static_cast<ScalarT>(this->outputValue(outputs, static_cast<Outputs>(n), RealT{0})));
        output_[n].initDerivative(ScalarT{0});
      }
      this->y_.setDataUpdated();
      this->yp_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
    {
      this->template parseInitialOutputs<Filter>(values);
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
    {
      return this->initializeOutputs(*this, values);
    }

    template <typename scalar_type, typename index_type>
    typename Filter<scalar_type, index_type>::Base::InitializationPortsT Filter<scalar_type, index_type>::initializationPorts()
    {
      typename Base::InitializationPortsT ports;
      for (const auto output : magic_enum::enum_values<Outputs>())
        if (output != Outputs::SIZE)
          ports.outputs.emplace(std::string(magic_enum::enum_name(output)), &outputSignal(output));
      return ports;
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::prepareInitialization(typename Base::InitialStateT& initial)
    {
      const auto outputs = this->template parseInitialOutputs<Filter>(initial.outputs(*this));
      for (size_t n = 0; n < output_.size(); ++n)
      {
        const auto value = this->outputValue(outputs, static_cast<Outputs>(n), RealT{0});
        initial.provide(output_[n], value);
        if (alias_[n])
          initial.provide(*alias_[n], value);
      }
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y, const ScalarT* yp, const ScalarT* voltage, const ScalarT*, ScalarT* f)
    {
      for (size_t n = 0; n < 3; ++n)
      {
        ScalarT source   = y[3 + n] - voltage[3 + n];
        ScalarT shunt    = y[6 + n] - y[n];
        ScalarT terminal = voltage[n] - y[3 + n];
        for (size_t k = 0; k < 3; ++k)
        {
          source   += Rs_[n][k] * y[k] + Ls_[n][k] * yp[k];
          shunt    += C_[n][k] * yp[3 + k];
          terminal += Rg_[n][k] * y[6 + k] + Lg_[n][k] * yp[6 + k];
        }
        f[n]     = source;
        f[3 + n] = shunt;
        f[6 + n] = terminal;
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();
      evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), this->y_ext_.data(), this->yp_ext_.data(), this->f_.getData());
      this->f_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::initializeMonitor()
    {
      using Mon = typename ModelDataT::MonitorableVariables;
      for (size_t n = 0; n < output_.size(); ++n)
        monitor_->set(static_cast<Mon>(n), [this, n]
                      { return output_[n].read(); });
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Filter<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }
  } // namespace EMT
} // namespace GridKit

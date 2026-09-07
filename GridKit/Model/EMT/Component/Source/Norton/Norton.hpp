#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class Norton : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using Base         = Component<ScalarT, IdxT>;
      using RealT        = typename Base::RealT;
      using SignalT      = Signal<ScalarT, IdxT>;
      using PhaseSignals = std::array<SignalT*, 3>;
      using YDataT       = VectorFitData<RealT, IdxT>;

      Norton(const YDataT& Y, PhaseSignals voltage, PhaseSignals incident, RealT scale)
        : voltage_(voltage), incident_(incident), admittance_(Y, scale)
      {
        if (Y.rows != 3 || Y.cols != 3 || Y.validate() != 0 || !std::isfinite(scale))
          throw std::invalid_argument("Norton: expected a finite three-phase admittance");
        this->equation_size_ = 3;
        this->size_          = 3 + admittance_.size();
        this->addOperator(&admittance_);
        for (auto& signal : shunt_)
          signal.claimProducer();
        admittance_.attachInput(voltage_[0], voltage_[1], voltage_[2]);
        admittance_.attachOutput(&shunt_[0], &shunt_[1], &shunt_[2]);
      }

      SignalT& voltage(size_t phase) const
      {
        return *voltage_.at(phase);
      }

      SignalT& outputSignal(size_t phase)
      {
        return shunt_.at(phase);
      }

      const SignalT& outputSignal(size_t phase) const
      {
        return shunt_.at(phase);
      }

      SignalT* inputSignal(size_t phase) const
      {
        return incident_.at(phase);
      }

      ScalarT incidentCurrent(size_t phase) const
      {
        const auto* signal = incident_.at(phase);
        return signal ? signal->read() : ScalarT{0};
      }

      ScalarT shuntCurrent(size_t phase) const
      {
        return shunt_.at(phase).read();
      }

      ScalarT current(size_t phase) const
      {
        return incidentCurrent(phase) - shunt_.at(phase).read();
      }

      int setGridKitComponentID(IdxT id) override
      {
        this->gridkit_component_id_ = id;
        return 0;
      }

      int allocate() override
      {
        if (!this->allocated_)
          this->allocateVectors(this->size_);
        this->tag_.resize(static_cast<size_t>(this->size_));
        this->variable_indices_.resize(static_cast<size_t>(this->size_));
        this->residual_indices_.resize(static_cast<size_t>(this->size_));
        for (IdxT p = 0; p < 3; ++p)
          this->bindSignal(shunt_[static_cast<size_t>(p)], p);
        const int status = this->allocateOperators();
        if (status != 0)
          return status;
        this->assignGlobalIndices(0);
        this->allocated_ = true;
        return 0;
      }

      int verify() const override
      {
        int errors = admittance_.verify();
        for (const auto* signal : incident_)
          if (signal && !signal->linked())
            ++errors;
        return errors;
      }

      int initializationOrder() const noexcept override
      {
        return 4;
      }

      int initializeState(const std::map<std::string, RealT>& values) override
      {
        if (!values.empty())
          throw std::invalid_argument("Norton shunt current is initialized from its admittance");
        return initialize();
      }

      int initialize()
      {
        const int status = admittance_.initialize();
        if (status != 0)
          return status;
        return initializeShunt();
      }

      int initializeShunt()
      {
        for (size_t p = 0; p < 3; ++p)
        {
          shunt_[p].init(admittance_.output(static_cast<IdxT>(p)));
          shunt_[p].initDerivative(ScalarT{0});
        }
        this->y_.setDataUpdated();
        this->yp_.setDataUpdated();
        return 0;
      }

      int initializeSteadyState(RealT omega)
      {
        std::array<RealT, 3> v, vp;
        for (size_t p = 0; p < 3; ++p)
        {
          v[p]  = static_cast<RealT>(voltage_[p]->read());
          vp[p] = static_cast<RealT>(voltage_[p]->readDerivative());
        }
        const int status = admittance_.initializeSteadyState(omega, v, vp);
        return status == 0 ? initializeShunt() : status;
      }

      int setAbsoluteTolerance(RealT tolerance) override
      {
        this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        return this->setAbsoluteToleranceOperators(tolerance);
      }

      int evaluateInternalResidual() override
      {
        auto* f = this->f_.getData();
        for (size_t p = 0; p < 3; ++p)
        {
          incident_value_[p] = incidentCurrent(p);
          f[p]               = -shunt_[p].read();
        }
        const int status = this->evaluateOperatorInternalResiduals();
        this->f_.setDataUpdated();
        return status;
      }

      int evaluateExternalResidual() override
      {
        for (size_t p = 0; p < 3; ++p)
          voltage_[p]->accumulateResidual(incident_value_[p] - shunt_[p].read());
        return this->evaluateOperatorExternalResiduals();
      }

      int evaluateResidual() override
      {
        const int status = evaluateInternalResidual();
        return status == 0 ? evaluateExternalResidual() : status;
      }

      int assembleJacobian(RealT y_scale, RealT yp_scale) override
      {
        std::array<typename SignalT::GradientT, 3> gradients;
        size_t                                     capacity = 6 + static_cast<size_t>(admittance_.jacobianCapacity()) * admittance_.externalJacobianExpansion();
        for (size_t p = 0; p < 3; ++p)
        {
          if (incident_[p])
            incident_[p]->appendGradient(gradients[p]);
          capacity += gradients[p].size();
        }
        if (capacity > jacobian_capacity_)
        {
          this->resetJacobianStructure();
          delete[] this->J_rows_buffer_;
          delete[] this->J_cols_buffer_;
          delete[] this->J_vals_buffer_;
          this->J_rows_buffer_ = new IdxT[capacity];
          this->J_cols_buffer_ = new IdxT[capacity];
          this->J_vals_buffer_ = new RealT[capacity];
          jacobian_capacity_   = capacity;
        }
        this->nnz_  = 0;
        auto append = [&](IdxT row, IdxT column, RealT value)
        {
          if (y_scale == ZERO<RealT>)
            return;
          value                   *= y_scale;
          const auto j             = this->nnz_++;
          this->J_rows_buffer_[j]  = row;
          this->J_cols_buffer_[j]  = column;
          this->J_vals_buffer_[j]  = value;
        };
        for (IdxT p = 0; p < 3; ++p)
        {
          append(this->getResidualIndex(p), this->getVariableIndex(p), -ONE<RealT>);
          append(voltage_[static_cast<size_t>(p)]->getResidualIndex(), this->getVariableIndex(p), -ONE<RealT>);
          for (const auto& [column, value] : gradients[static_cast<size_t>(p)])
            append(voltage_[static_cast<size_t>(p)]->getResidualIndex(), column, value);
        }
        const int status = this->evaluateOperatorJacobians(y_scale, yp_scale);
        if (status != 0)
          return status;
        this->appendOperatorJacobians();
        return this->constructCoo();
      }

    private:
      PhaseSignals             voltage_, incident_;
      std::array<ScalarT, 3>   incident_value_{};
      std::array<SignalT, 3>   shunt_;
      VectorFit<ScalarT, IdxT> admittance_;
      size_t                   jacobian_capacity_{0};
    };
  } // namespace EMT
} // namespace GridKit

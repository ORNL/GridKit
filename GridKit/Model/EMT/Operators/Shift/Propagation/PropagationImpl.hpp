#pragma once

#include <GridKit/Model/EMT/Operators/Shift/Propagation/Propagation.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Propagation<scalar_type, index_type>::Propagation(const ModelDataT& data)
      : data_(data), filtered_(data.modes.size() * static_cast<size_t>(data.K)), output_(static_cast<size_t>(data.K)), delay_(delayData(data))
    {
      if (data.validate())
        throw std::invalid_argument("Propagation: expected stable proper square modes with positive delays");
      const auto K         = static_cast<size_t>(data.K);
      this->equation_size_ = 0;
      for (size_t m = 0; m < data.modes.size(); ++m)
      {
        auto& fit = *fits_.emplace_back(std::make_unique<VectorFitT>(data.modes[m].H));
        this->addOperator(&fit);
        this->size_ += fit.size();
        for (size_t k = 0; k < K; ++k)
          filtered_[m * K + k].setComputed(
              [this, m, k]
              { return fits_[m]->output(static_cast<IdxT>(k)); },
              [this, m, k](typename SignalT::GradientT& gradient, RealT scale)
              { fits_[m]->appendOutputGradient(static_cast<IdxT>(k), gradient, scale); });
      }
      std::vector<SignalT*> signals;
      for (auto& signal : filtered_)
        signals.push_back(&signal);
      delay_.attachInput(signals);
      delay_.setInputDerivative([this, K](size_t channel)
                                { return static_cast<RealT>(fits_[channel / K]->outputDerivative(static_cast<IdxT>(channel % K))); });
      this->addOperator(&delay_);
      this->size_ += delay_.size();
      for (size_t k = 0; k < K; ++k)
      {
        output_[k].claimProducer();
        output_[k].setComputed(
            [this, k, K]
            {
              ScalarT value{0};
              for (size_t m = 0; m < fits_.size(); ++m)
                value += delay_.outputSignal(static_cast<IdxT>(m * K + k)).read();
              return value;
            },
            [this, k, K](typename SignalT::GradientT& gradient, RealT scale)
            {
              for (size_t m = 0; m < fits_.size(); ++m)
                delay_.outputSignal(static_cast<IdxT>(m * K + k)).appendGradient(gradient, scale);
            });
      }
    }

    template <typename scalar_type, typename index_type>
    int Propagation<scalar_type, index_type>::allocate()
    {
      if (!this->allocated_)
        this->allocateVectors(this->size_);
      this->tag_.resize(static_cast<size_t>(this->size_));
      this->variable_indices_.resize(static_cast<size_t>(this->size_));
      this->residual_indices_.resize(static_cast<size_t>(this->size_));
      const int status = this->allocateOperators();
      if (status != 0)
        return status;
      this->assignGlobalIndices(0);
      this->allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Propagation<scalar_type, index_type>::verify() const
    {
      int errors = data_.validate() + delay_.verify();
      for (const auto& fit : fits_)
        errors += fit->verify();
      return errors;
    }

    template <typename scalar_type, typename index_type>
    int Propagation<scalar_type, index_type>::initializeSteadyState(RealT omega, std::span<const RealT> u, std::span<const RealT> udot, RealT time)
    {
      if (u.size() != output_.size() || udot.size() != u.size())
        throw std::invalid_argument("Propagation: prehistory dimension mismatch");
      std::vector<RealT> values(filtered_.size()), slopes(filtered_.size());
      for (size_t m = 0; m < fits_.size(); ++m)
      {
        const int status = fits_[m]->initializeSteadyState(omega, u, udot);
        if (status != 0)
          return status;
        MatrixT re, im;
        fits_[m]->transfer(omega, re, im);
        for (size_t n = 0; n < u.size(); ++n)
          for (size_t k = 0; k < u.size(); ++k)
          {
            values[m * u.size() + n] += re[n][k] * u[k] + (omega == RealT{0} ? RealT{0} : im[n][k] * udot[k] / omega);
            slopes[m * u.size() + n] += re[n][k] * udot[k] - omega * im[n][k] * u[k];
          }
      }
      return delay_.initializeSteadyState(omega, values, slopes, time);
    }

    template <typename scalar_type, typename index_type>
    int Propagation<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      return this->setAbsoluteToleranceOperators(tolerance);
    }

    template <typename scalar_type, typename index_type>
    int Propagation<scalar_type, index_type>::evaluateInternalResidual()
    {
      return this->evaluateOperatorInternalResiduals();
    }

    template <typename scalar_type, typename index_type>
    int Propagation<scalar_type, index_type>::evaluateResidual()
    {
      const int status = evaluateInternalResidual();
      return status == 0 ? this->evaluateExternalResidual() : status;
    }

    template <typename scalar_type, typename index_type>
    int Propagation<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      const int status = this->evaluateOperatorJacobians(y_scale, yp_scale);
      if (status != 0)
        return status;
      size_t capacity = 0;
      for (auto* op : this->operators_)
        if (op->getCooJacobian())
          capacity += static_cast<size_t>(op->getCooJacobian()->getNnz());
      // Growing parent storage must not invalidate the child triplets just evaluated.
      if (capacity > capacity_)
      {
        delete this->coo_jac_;
        this->coo_jac_ = nullptr;
        delete[] this->J_rows_buffer_;
        delete[] this->J_cols_buffer_;
        delete[] this->J_vals_buffer_;
        this->J_rows_buffer_ = new IdxT[capacity];
        this->J_cols_buffer_ = new IdxT[capacity];
        this->J_vals_buffer_ = new RealT[capacity];
        capacity_            = capacity;
      }
      this->nnz_ = 0;
      this->appendOperatorJacobians();
      return this->constructCoo();
    }

    template <typename scalar_type, typename index_type>
    void Propagation<scalar_type, index_type>::attachInput(const std::vector<SignalT*>& input)
    {
      if (this->allocated_ || input.size() != output_.size())
        throw std::invalid_argument("Propagation: invalid input binding");
      for (auto& fit : fits_)
        fit->attachInput(input);
    }

    template <typename scalar_type, typename index_type>
    typename Propagation<scalar_type, index_type>::SignalT& Propagation<scalar_type, index_type>::outputSignal(IdxT channel)
    {
      return output_.at(static_cast<size_t>(channel));
    }

    template <typename scalar_type, typename index_type>
    void Propagation<scalar_type, index_type>::transfer(RealT omega, MatrixT& re, MatrixT& im) const
    {
      re = MatrixT(output_.size(), output_.size());
      im = MatrixT(output_.size(), output_.size());
      for (size_t m = 0; m < fits_.size(); ++m)
      {
        MatrixT r, i;
        fits_[m]->transfer(omega, r, i);
        const RealT c = std::cos(omega * data_.modes[m].tau), s = std::sin(omega * data_.modes[m].tau);
        for (size_t n = 0; n < output_.size(); ++n)
          for (size_t k = 0; k < output_.size(); ++k)
          {
            re[n][k] += c * r[n][k] + s * i[n][k];
            im[n][k] += c * i[n][k] - s * r[n][k];
          }
      }
    }

    template <typename scalar_type, typename index_type>
    DelayData<typename Propagation<scalar_type, index_type>::RealT, index_type>
    Propagation<scalar_type, index_type>::delayData(const ModelDataT& data)
    {
      if (data.validate())
        throw std::invalid_argument("Propagation: invalid mode coefficients");
      DelayData<RealT, IdxT> delay;
      delay.M          = static_cast<IdxT>(data.modes.size()) * data.K;
      delay.limit_step = data.limit_step;
      for (const auto& mode : data.modes)
        for (IdxT k = 0; k < data.K; ++k)
          delay.tau.push_back(mode.tau);
      return delay;
    }
  } // namespace EMT
} // namespace GridKit

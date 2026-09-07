#pragma once

#include <GridKit/Model/EMT/Operators/Shift/Delay/Delay.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Delay<scalar_type, index_type>::Delay(const ModelDataT& data)
      : data_(data), input_(data.tau.size()), output_(data.tau.size()), constant_(data.tau.size()), coefficient_(data.tau.size()), gradients_(data.tau.size())
    {
      if (data.validate())
        throw std::invalid_argument("Delay: expected positive finite channel delays");
      this->equation_size_ = this->size_ = data.M;
      tau_max_                           = *std::max_element(data.tau.begin(), data.tau.end());
      for (auto& signal : output_)
        signal.claimProducer();
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::attachInput(const std::vector<SignalT*>& input)
    {
      if (this->allocated_ || input.size() != input_.size())
        throw std::invalid_argument("Delay: invalid input binding");
      input_ = input;
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::setInputDerivative(std::function<RealT(size_t)> derivative)
    {
      derivative_ = std::move(derivative);
    }

    template <typename scalar_type, typename index_type>
    typename Delay<scalar_type, index_type>::SignalT& Delay<scalar_type, index_type>::outputSignal(IdxT channel)
    {
      return output_.at(static_cast<size_t>(channel));
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::setPrehistory(RealT time, Prehistory prehistory)
    {
      if (!std::isfinite(time) || !prehistory)
        throw std::invalid_argument("Delay: a finite initial time and prehistory are required");
      origin_     = time;
      prehistory_ = std::move(prehistory);
      resetHistory();
    }

    template <typename scalar_type, typename index_type>
    int Delay<scalar_type, index_type>::initializeSteadyState(RealT omega, std::span<const RealT> u, std::span<const RealT> udot, RealT time)
    {
      if (!std::isfinite(omega) || omega < RealT{0} || u.size() != input_.size() || udot.size() != u.size())
        throw std::invalid_argument("Delay: invalid harmonic prehistory");
      std::vector<RealT> value(u.begin(), u.end()), slope(udot.begin(), udot.end());
      for (size_t k = 0; k < u.size(); ++k)
        if (!std::isfinite(u[k]) || !std::isfinite(udot[k]) || (omega == RealT{0} && udot[k] != RealT{0}))
          throw std::invalid_argument("Delay: invalid harmonic prehistory value");
      setPrehistory(time, [value, slope, omega, time](size_t k, RealT t)
                    {
                      if (omega == RealT{0})
                        return std::pair{value[k], RealT{0}};
                      const RealT angle = omega * (t - time), c = std::cos(angle), s = std::sin(angle);
                      return std::pair{value[k] * c + slope[k] * s / omega,
                                       -omega * value[k] * s + slope[k] * c}; });
      return initialize();
    }

    template <typename scalar_type, typename index_type>
    int Delay<scalar_type, index_type>::allocate()
    {
      if (!this->allocated_)
        this->allocateVectors(this->size_);
      this->tag_.resize(input_.size());
      this->variable_indices_.resize(input_.size());
      this->residual_indices_.resize(input_.size());
      this->assignGlobalIndices(0);
      this->allocateExternalVectors(data_.M, 0);
      for (IdxT k = 0; k < data_.M; ++k)
      {
        this->bindSignal(output_[static_cast<size_t>(k)], k);
        this->setExternalVariableSignal(k, input_[static_cast<size_t>(k)]);
      }
      this->allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Delay<scalar_type, index_type>::verify() const
    {
      int errors = data_.validate();
      for (auto* signal : input_)
        errors += !signal || (this->allocated_ && (!signal->linked() || (!derivative_ && !signal->derivativeLinked())));
      return errors;
    }

    template <typename scalar_type, typename index_type>
    int Delay<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::transfer(RealT omega, MatrixT& re, MatrixT& im) const
    {
      if (!std::isfinite(omega))
        throw std::invalid_argument("Delay: frequency must be finite");
      re = MatrixT(input_.size(), input_.size());
      im = MatrixT(input_.size(), input_.size());
      for (size_t k = 0; k < input_.size(); ++k)
      {
        re[k][k] = std::cos(omega * data_.tau[k]);
        im[k][k] = -std::sin(omega * data_.tau[k]);
      }
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::acceptStep(RealT time)
    {
      if (!prehistory_ || knots_.empty() || time < knots_.back().time)
        throw std::logic_error("Delay: accepted history must advance from its prehistory");
      const bool restart = time == knots_.back().time;
      Knot       knot{time, {}, {}, !restart};
      knot.value.reserve(input_.size());
      knot.slope.reserve(input_.size());
      for (size_t k = 0; k < input_.size(); ++k)
      {
        const RealT value = static_cast<RealT>(input_[k]->read());
        const RealT slope = derivative_ ? derivative_(k) : static_cast<RealT>(input_[k]->readDerivative());
        if (!std::isfinite(value) || !std::isfinite(slope))
          throw std::runtime_error("Delay: nonfinite accepted input");
        if (restart && std::abs(value - knots_.back().value[k]) > RealT{1e-10} * (ONE<RealT> + std::abs(value)))
          discontinuities_.insert(time + data_.tau[k]);
        knot.value.push_back(value);
        knot.slope.push_back(slope);
      }
      // Keep both sides of an event. IDACalcIC does not refresh algebraic slopes.
      knots_.push_back(std::move(knot));
      const RealT oldest = time - tau_max_;
      while (knots_.size() > 3 && knots_[1].time < oldest)
        knots_.pop_front();
    }

    template <typename scalar_type, typename index_type>
    typename Delay<scalar_type, index_type>::RealT Delay<scalar_type, index_type>::maximumStepSize() const
    {
      return data_.limit_step ? *std::min_element(data_.tau.begin(), data_.tau.end())
                              : std::numeric_limits<RealT>::infinity();
    }

    template <typename scalar_type, typename index_type>
    int Delay<scalar_type, index_type>::evaluateInternalResidual()
    {
      if (knots_.empty())
        throw std::logic_error("Delay: explicit prehistory is required before residual evaluation");
      for (size_t k = 0; k < input_.size(); ++k)
        this->f_.getData()[k] = -output_[k].read() + constant_[k] + coefficient_[k] * input_[k]->read();
      this->f_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Delay<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    int Delay<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT)
    {
      size_t capacity = input_.size();
      for (size_t k = 0; k < input_.size(); ++k)
      {
        gradients_[k].clear();
        if (y_scale != RealT{0})
          input_[k]->appendGradient(gradients_[k]);
        capacity += gradients_[k].size();
      }
      this->reserveJacobian(capacity);
      this->nnz_ = 0;
      if (y_scale != RealT{0})
        for (IdxT k = 0; k < data_.M; ++k)
        {
          const auto row = this->getResidualIndex(k);
          this->appendJacobian(row, this->getVariableIndex(k), -y_scale);
          for (const auto& [column, value] : gradients_[static_cast<size_t>(k)])
            this->appendJacobian(row, column, y_scale * coefficient_[static_cast<size_t>(k)] * value);
        }
      return this->constructCoo();
    }

    template <typename scalar_type, typename index_type>
    int Delay<scalar_type, index_type>::initialize()
    {
      if (!this->allocated_ || !prehistory_)
        throw std::logic_error("Delay: allocation and explicit prehistory are required");
      for (size_t k = 0; k < input_.size(); ++k)
      {
        const auto [value, slope] = prehistory_(k, origin_ - data_.tau[k]);
        output_[k].init(static_cast<ScalarT>(value));
        output_[k].initDerivative(static_cast<ScalarT>(slope));
      }
      this->y_.setDataUpdated();
      this->yp_.setDataUpdated();
      updateTime(origin_, ONE<RealT>);
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::resetHistory()
    {
      knots_.clear();
      discontinuities_.clear();
      right_limit_time_ = origin_;
      if (prehistory_)
      {
        Knot knot{origin_, {}, {}, true};
        for (size_t k = 0; k < input_.size(); ++k)
        {
          const auto [value, slope] = prehistory_(k, origin_);
          if (!std::isfinite(value) || !std::isfinite(slope))
            throw std::invalid_argument("Delay: nonfinite prehistory");
          knot.value.push_back(value);
          knot.slope.push_back(slope);
        }
        knots_.push_back(std::move(knot));
      }
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::updateTime(RealT time, RealT alpha)
    {
      Base::updateTime(time, alpha);
      if (knots_.empty())
        return;
      const auto& last = knots_.back();
      for (size_t k = 0; k < input_.size(); ++k)
      {
        const RealT xi  = time - data_.tau[k];
        coefficient_[k] = ZERO<RealT>;
        if (xi <= last.time)
          constant_[k] = history(k, xi);
        else
        {
          const RealT h = time - last.time, theta = (xi - last.time) / h;
          coefficient_[k] = last.smooth ? theta * theta : theta;
          constant_[k]    = (ONE<RealT> - coefficient_[k]) * last.value[k];
          if (last.smooth)
            constant_[k] += h * theta * (ONE<RealT> - theta) * last.slope[k];
        }
      }
    }

    template <typename scalar_type, typename index_type>
    typename Delay<scalar_type, index_type>::RealT Delay<scalar_type, index_type>::history(size_t k, RealT time) const
    {
      if (time < origin_)
        return prehistory_(k, time).first;
      // Select the right limit at a duplicate event knot.
      const RealT epsilon   = RealT{32} * std::numeric_limits<RealT>::epsilon() * std::max(std::abs(this->time_), tau_max_);
      const auto  candidate = std::lower_bound(knots_.begin(), knots_.end(), time - epsilon, [](const Knot& knot, RealT t)
                                              { return knot.time < t; });
      if (candidate != knots_.end() && std::abs(candidate->time - time) <= epsilon)
      {
        time = candidate->time;
        if (this->time_ > right_limit_time_ + epsilon)
          return candidate->value[k];
      }
      const auto upper = std::upper_bound(knots_.begin(), knots_.end(), time, [](RealT t, const Knot& knot)
                                          { return t < knot.time; });
      if (upper == knots_.begin())
        throw std::logic_error("Delay: history lookup predates retained knots");
      const auto& left = *std::prev(upper);
      if (upper == knots_.end() || time == left.time)
        return left.value[k];
      const auto& right = *upper;
      const RealT h = right.time - left.time, theta = (time - left.time) / h;
      if (!left.smooth)
        return (ONE<RealT> - theta) * left.value[k] + theta * right.value[k];
      const RealT a = ONE<RealT> - theta;
      return a * a * (ONE<RealT> + TWO<RealT> * theta) * left.value[k]
             + theta * a * a * h * left.slope[k]
             + theta * theta * (THREE<RealT> - TWO<RealT> * theta) * right.value[k]
             - theta * theta * a * h * right.slope[k];
    }

    template <typename scalar_type, typename index_type>
    typename Delay<scalar_type, index_type>::RealT Delay<scalar_type, index_type>::nextDiscontinuityTime(RealT after) const
    {
      const RealT epsilon = RealT{32} * std::numeric_limits<RealT>::epsilon() * std::max(std::abs(after), tau_max_);
      const auto  next    = discontinuities_.upper_bound(after + epsilon);
      return next == discontinuities_.end() ? std::numeric_limits<RealT>::infinity() : *next;
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::beginDiscontinuity(RealT time)
    {
      right_limit_time_   = time;
      const RealT epsilon = RealT{32} * std::numeric_limits<RealT>::epsilon() * std::max(std::abs(time), tau_max_);
      discontinuities_.erase(discontinuities_.begin(), discontinuities_.upper_bound(time + epsilon));
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::reserveJacobian(size_t capacity)
    {
      capacity = std::max(size_t{1}, capacity);
      if (capacity <= jacobian_capacity_)
        return;
      this->resetJacobianStructure();
      delete[] this->J_rows_buffer_;
      delete[] this->J_cols_buffer_;
      delete[] this->J_vals_buffer_;
      this->J_rows_buffer_ = new IdxT[capacity];
      this->J_cols_buffer_ = new IdxT[capacity];
      this->J_vals_buffer_ = new RealT[capacity];
      jacobian_capacity_   = capacity;
    }

    template <typename scalar_type, typename index_type>
    void Delay<scalar_type, index_type>::appendJacobian(IdxT row, IdxT column, RealT value)
    {
      const auto j            = this->nnz_++;
      this->J_rows_buffer_[j] = row;
      this->J_cols_buffer_[j] = column;
      this->J_vals_buffer_[j] = value;
    }

  } // namespace EMT
} // namespace GridKit

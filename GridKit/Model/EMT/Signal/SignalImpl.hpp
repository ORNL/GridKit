/**
 * @file Signal model implementation.
 */
#include <stdexcept>
#include <utility>

#include <GridKit/Model/EMT/Signal/Signal.hpp>
#include <GridKit/Model/EMT/Signal/SignalData.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Signal<scalar_type, index_type>::Signal()
    {
    }

    template <typename scalar_type, typename index_type>
    Signal<scalar_type, index_type>::Signal(std::string id)
      : id_(std::move(id))
    {
    }

    template <typename scalar_type, typename index_type>
    Signal<scalar_type, index_type>::Signal(const SignalData<RealT, IdxT>& data)
      : Signal(data.id)
    {
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::set(ScalarT* signal, IdxT* variable_index)
    {
      set(signal, nullptr, nullptr, variable_index, nullptr);
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::set(ScalarT* signal,
                                              ScalarT* derivative,
                                              ScalarT* residual,
                                              IdxT*    variable_index,
                                              IdxT*    residual_index)
    {
      producer_claimed_ = true;
      signal_           = signal;
      derivative_       = derivative;
      residual_         = residual;
      variable_index_   = variable_index;
      residual_index_   = residual_index;
      value_            = {};
      gradient_         = {};
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::setComputed(
        std::function<ScalarT()>               value,
        std::function<void(GradientT&, RealT)> gradient)
    {
      if (!value || !gradient)
      {
        throw std::invalid_argument("A computed signal requires a value and gradient");
      }
      set(nullptr, nullptr);
      value_    = std::move(value);
      gradient_ = std::move(gradient);
    }

    template <typename scalar_type, typename index_type>
    bool Signal<scalar_type, index_type>::computed() const
    {
      return static_cast<bool>(value_);
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::appendGradient(GradientT& gradient, RealT scale) const
    {
      if (computed())
      {
        if (evaluating_)
        {
          throw std::logic_error("Cyclic computed signal \"" + id_ + "\"");
        }
        evaluating_ = true;
        try
        {
          gradient_(gradient, scale);
        }
        catch (...)
        {
          evaluating_ = false;
          throw;
        }
        evaluating_ = false;
      }
      else if (linked() && getVariableIndex() != INVALID_INDEX<IdxT>)
      {
        gradient.emplace_back(getVariableIndex(), scale);
      }
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::claimProducer()
    {
      if (producer_claimed_)
      {
        throw std::logic_error("Signal \"" + id_ + "\" has more than one producer");
      }
      producer_claimed_ = true;
    }

    template <typename scalar_type, typename index_type>
    bool Signal<scalar_type, index_type>::hasProducer() const
    {
      return producer_claimed_;
    }

    template <typename scalar_type, typename index_type>
    bool Signal<scalar_type, index_type>::linked() const
    {
      return computed() || ((signal_) && (variable_index_));
    }

    template <typename scalar_type, typename index_type>
    bool Signal<scalar_type, index_type>::derivativeLinked() const
    {
      return derivative_ != nullptr;
    }

    template <typename scalar_type, typename index_type>
    bool Signal<scalar_type, index_type>::residualLinked() const
    {
      return (residual_) && (residual_index_);
    }

    template <typename scalar_type, typename index_type>
    scalar_type Signal<scalar_type, index_type>::read() const
    {
      if (computed())
      {
        if (evaluating_)
        {
          throw std::logic_error("Cyclic computed signal \"" + id_ + "\"");
        }
        evaluating_ = true;
        try
        {
          auto value  = value_();
          evaluating_ = false;
          return value;
        }
        catch (...)
        {
          evaluating_ = false;
          throw;
        }
      }
      if (!linked())
      {
        throw std::logic_error("Unlinked signal \"" + id_ + "\"");
      }
      return *signal_;
    }

    template <typename scalar_type, typename index_type>
    scalar_type Signal<scalar_type, index_type>::readDerivative() const
    {
      if (!derivativeLinked())
      {
        throw std::logic_error("Signal \"" + id_ + "\" has no derivative");
      }
      return *derivative_;
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::accumulateResidual(ScalarT value)
    {
      *residual_ += value;
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::init(ScalarT signal)
    {
      if (computed())
      {
        throw std::logic_error("Cannot initialize a computed signal");
      }
      *signal_ = signal;
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::initDerivative(ScalarT derivative)
    {
      if (!derivativeLinked())
      {
        throw std::logic_error("Signal \"" + id_ + "\" has no derivative");
      }
      *derivative_ = derivative;
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::markDerivativeCoupling()
    {
      if (computed())
      {
        throw std::logic_error("Computed algebraic signals do not expose derivatives");
      }
      derivative_coupling_ = true;
    }

    template <typename scalar_type, typename index_type>
    bool Signal<scalar_type, index_type>::hasDerivativeCoupling() const
    {
      return derivative_coupling_;
    }
  } // namespace EMT
} // namespace GridKit

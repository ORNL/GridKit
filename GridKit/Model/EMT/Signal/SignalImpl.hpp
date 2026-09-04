/**
 * @file Signal model implementation.
 */
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
      signal_         = signal;
      variable_index_ = variable_index;
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::set(ScalarT* signal,
                                              ScalarT* derivative,
                                              ScalarT* residual,
                                              IdxT*    variable_index,
                                              IdxT*    residual_index)
    {
      signal_         = signal;
      derivative_     = derivative;
      residual_       = residual;
      variable_index_ = variable_index;
      residual_index_ = residual_index;
    }

    template <typename scalar_type, typename index_type>
    bool Signal<scalar_type, index_type>::linked() const
    {
      return (signal_) && (variable_index_);
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
      return *signal_;
    }

    template <typename scalar_type, typename index_type>
    scalar_type Signal<scalar_type, index_type>::readDerivative() const
    {
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
      *signal_ = signal;
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::initDerivative(ScalarT derivative)
    {
      *derivative_ = derivative;
    }

    template <typename scalar_type, typename index_type>
    void Signal<scalar_type, index_type>::markDerivativeCoupling()
    {
      derivative_coupling_ = true;
    }

    template <typename scalar_type, typename index_type>
    bool Signal<scalar_type, index_type>::hasDerivativeCoupling() const
    {
      return derivative_coupling_;
    }
  } // namespace EMT
} // namespace GridKit

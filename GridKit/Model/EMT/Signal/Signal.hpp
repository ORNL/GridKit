#pragma once

#include <array>
#include <string>

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct SignalData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /*!
     * @brief Signal model implementation base class.
     *
     * A fully bound signal exposes the value, derivative, and residual row of
     * one variable owned by another component, along with the global variable
     * and residual indices. Consumers read the value and derivative through
     * the signal and accumulate external residual contributions into the
     * owner's residual row.
     */
    template <typename scalar_type, typename index_type>
    class Signal
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      Signal();
      explicit Signal(std::string id);
      explicit Signal(const SignalData<RealT, IdxT>& data);

      virtual ~Signal() = default;

      void set(ScalarT* signal, IdxT* variable_index);
      void set(ScalarT* signal,
               ScalarT* derivative,
               ScalarT* residual,
               IdxT*    variable_index,
               IdxT*    residual_index);

      bool linked() const;
      bool derivativeLinked() const;
      bool residualLinked() const;

      ScalarT read() const;
      ScalarT readDerivative() const;
      void    accumulateResidual(ScalarT value);
      void    init(ScalarT signal_in);
      void    initDerivative(ScalarT derivative_in);

      void markDerivativeCoupling();
      bool hasDerivativeCoupling() const;

      const std::string& id() const noexcept
      {
        return id_;
      }

      IdxT getVariableIndex() const
      {
        return *variable_index_;
      }

      IdxT getResidualIndex() const
      {
        return *residual_index_;
      }

    private:
      ScalarT*    signal_{nullptr};
      ScalarT*    derivative_{nullptr};
      ScalarT*    residual_{nullptr};
      std::string id_;

    protected:
      IdxT* variable_index_{nullptr};
      IdxT* residual_index_{nullptr};
      bool  derivative_coupling_{false};
    };

    /*!
     * @brief Three-phase electrical port bundling one signal per phase.
     *
     * Owned by value by the component that owns the phase variables and
     * residual rows.
     */
    template <typename scalar_type, typename index_type>
    struct Port3
    {
      using SignalT = Signal<scalar_type, index_type>;

      std::array<SignalT, 3> signals{};

      SignalT* a()
      {
        return &signals[0];
      }

      SignalT* b()
      {
        return &signals[1];
      }

      SignalT* c()
      {
        return &signals[2];
      }
    };

  } // namespace EMT
} // namespace GridKit

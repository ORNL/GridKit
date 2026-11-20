/**
 * @file SignalImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT signal.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Signal/Signal.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct `N` signals.
     *
     * Layout: internal variable `s` of size `N`, with no terminal, input or
     * output.
     */
    template <typename scalar_type, typename index_type>
    Signal<scalar_type, index_type>::Signal(const ModelDataT& data)
      : N_(data.N),
        held_(data.held)
    {
      this->setLayout({N_}, {}, {}, {});
    }

    /**
     * @brief Check that `N` is positive.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Signal<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Signal");
      checks.check(N_ > 0, "N must be positive");
      return checks;
    }

    /**
     * @brief Write the residual of `s`: minus its derivative if it is held,
     *        and otherwise minus `s`, to which the output it follows is added.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void Signal<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T*, const T*, T* f) const
    {
      const IdxT S = this->local(SignalInternalVariables::s);

      const T* s = y;
      if (held_)
      {
        s = yp;
      }
      for (IdxT n = 0; n < N_; ++n)
      {
        f[S + n] = -s[S + n];
      }
    }
  } // namespace EMT
} // namespace GridKit

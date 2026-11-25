/**
 * @file BusImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT bus.
 *
 * The model equations are in the `README.md` of this directory.
 */

#pragma once

#include <algorithm>

#include <GridKit/Model/EMT/Bus/Bus.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a bus with `N` phases.
     *
     * Layout: internal variable `v` of size `N`, with no terminal, input or
     * output. The current balance holds only what the devices add at their
     * terminals.
     */
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : N_(data.N),
        kv_(data.kv)
    {
      this->setLayout({N_}, {}, {}, {});
    }

    /**
     * @brief Check that `N` and `kv` are positive.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Bus<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Bus");
      checks.check(N_ > 0, "N must be positive");
      checks.check(kv_ > ZERO<RealT>, "kv must be positive");
      return checks;
    }

    /**
     * @brief Set the tolerances of the voltages to `rel_tol` times the bus
     *        base voltage.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre The system has bound this bus.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      std::ranges::fill(this->absoluteTolerance(BusInternalVariables::v), rel_tol * KILO<RealT> * kv_);
      return 0;
    }
  } // namespace EMT
} // namespace GridKit

/**
 * @file BusImpl.hpp
 * @brief Definition of an optimal power flow bus.
 */

#pragma once

#include <format>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModelImpl.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Construct a bus from data
     *
     * The power balance is an equality, the squared voltage magnitude has the
     * voltage limits, and a reference bus limits its imaginary voltage to zero.
     */
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : ComponentModelT(std::format("bus_{}", data.number), {}),
        reference_(data.reference),
        Vmin_(parameter(data, BusParameters::Vmin, ZERO<RealT>)),
        Vmax_(parameter(data, BusParameters::Vmax, UNBOUNDED<RealT>))
    {
      constraint_lower_[0] = ZERO<RealT>;
      constraint_upper_[0] = ZERO<RealT>;
      constraint_lower_[1] = ZERO<RealT>;
      constraint_upper_[1] = ZERO<RealT>;
      constraint_upper_[2] = Vmax_ * Vmax_;
      if (Vmin_ > ZERO<RealT>)
      {
        constraint_lower_[2] = Vmin_ * Vmin_;
      }

      if (reference_)
      {
        variable_lower_[1] = ZERO<RealT>;
        variable_upper_[1] = ZERO<RealT>;
      }
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Bus<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Bus");
      checks.check(Vmin_ >= ZERO<RealT>, "Vmin must be nonnegative");
      checks.check(Vmax_ >= Vmin_, "Vmax must not be less than Vmin");
      return checks;
    }

    /**
     * @brief Flat start
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize(ScalarT* x)
    {
      x[variable_indices_[0]] = ONE<RealT>;
      x[variable_indices_[1]] = ZERO<RealT>;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline typename Bus<scalar_type, index_type>::ScalarT
    Bus<scalar_type, index_type>::objective([[maybe_unused]] const ScalarT* x) const
    {
      return ScalarT{ZERO<RealT>};
    }

    /**
     * @brief Squared voltage magnitude
     *
     * Devices add their injections to the balance rows, which the bus leaves
     * at zero.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void
    Bus<scalar_type, index_type>::constraints(const ScalarT* x, ScalarT* g) const
    {
      const ScalarT vr = x[0];
      const ScalarT vi = x[1];

      g[2] = vr * vr + vi * vi;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit

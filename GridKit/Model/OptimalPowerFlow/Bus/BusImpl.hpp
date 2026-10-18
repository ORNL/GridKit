/**
 * @file BusImpl.hpp
 * @brief Definition of an optimal power flow bus.
 */

#pragma once

#include <cmath>
#include <format>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : ComponentT(std::type_identity<Bus>{}, std::format("bus {}", data.number)),
        reference_(data.reference),
        Vmin_(parameter(data, BusParameters::Vmin, ZERO<RealT>)),
        Vmax_(parameter(data, BusParameters::Vmax, UNBOUNDED<RealT>)),
        G_(parameter(data, BusParameters::G, ZERO<RealT>)),
        B_(parameter(data, BusParameters::B, ZERO<RealT>))
    {
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Bus<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Bus");
      checks.check(Vmin_ >= ZERO<RealT>, "Vmin must be nonnegative");
      checks.check(Vmax_ >= Vmin_, "Vmax must not be less than Vmin");
      checks.check(std::isfinite(G_), "G must be finite");
      checks.check(std::isfinite(B_), "B must be finite");
      return checks;
    }

    /**
     * @brief The voltage limits bound the squared magnitude, and a reference
     * bus fixes its imaginary voltage at zero
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::allocate()
    {
      d_upper_.getData()[0] = Vmax_ * Vmax_;
      if (Vmin_ > ZERO<RealT>)
      {
        d_lower_.getData()[0] = Vmin_ * Vmin_;
      }

      if (reference_)
      {
        x_lower_.getData()[1] = ZERO<RealT>;
        x_upper_.getData()[1] = ZERO<RealT>;
      }
      return 0;
    }

    /**
     * @brief Flat start
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize()
    {
      ScalarT* x = x_.getData();
      x[0]       = ONE<RealT>;
      x[1]       = ZERO<RealT>;
      return 0;
    }

    /**
     * @brief Active and reactive power into the bus from its shunt; the
     * ports add the rest of the balance
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Bus<scalar_type, index_type>::internalEqualities(const ScalarT* x, ScalarT* c) const
    {
      const ScalarT v2 = x[0] * x[0] + x[1] * x[1];

      c[0] = -G_ * v2;
      c[1] = B_ * v2;
    }

    /**
     * @brief Squared voltage magnitude
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Bus<scalar_type, index_type>::internalInequalities(const ScalarT* x, ScalarT* d) const
    {
      const ScalarT vr = x[0];
      const ScalarT vi = x[1];

      d[0] = vr * vr + vi * vi;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit

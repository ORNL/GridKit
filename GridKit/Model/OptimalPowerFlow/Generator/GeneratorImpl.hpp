/**
 * @file GeneratorImpl.hpp
 * @brief Definition of an optimal power flow generator.
 */

#pragma once

#include <cmath>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Generator<scalar_type, index_type>::Generator(const ModelDataT& data)
      : ComponentT(std::type_identity<Generator>{}, data.id),
        in_service_(data.in_service),
        Pmin_(parameter(data, GeneratorParameters::Pmin, -UNBOUNDED<RealT>)),
        Pmax_(parameter(data, GeneratorParameters::Pmax, UNBOUNDED<RealT>)),
        Qmin_(parameter(data, GeneratorParameters::Qmin, -UNBOUNDED<RealT>)),
        Qmax_(parameter(data, GeneratorParameters::Qmax, UNBOUNDED<RealT>)),
        c0_(parameter(data, GeneratorParameters::c0, ZERO<RealT>)),
        c1_(parameter(data, GeneratorParameters::c1, ZERO<RealT>)),
        c2_(parameter(data, GeneratorParameters::c2, ZERO<RealT>))
    {
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Generator<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Generator");
      checks.check(Pmin_ <= Pmax_, "Pmin must not exceed Pmax");
      checks.check(Qmin_ <= Qmax_, "Qmin must not exceed Qmax");
      checks.check(std::isfinite(c0_), "c0 must be finite");
      checks.check(std::isfinite(c1_), "c1 must be finite");
      checks.check(std::isfinite(c2_), "c2 must be finite");
      return checks;
    }

    /**
     * @brief The power limits bound the power, which is zero out of service
     */
    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::allocate()
    {
      const bool in_service = in_service_ > ZERO<RealT>;
      RealT*     lower      = x_lower_.getData();
      RealT*     upper      = x_upper_.getData();
      lower[0]              = in_service ? Pmin_ : ZERO<RealT>;
      lower[1]              = in_service ? Qmin_ : ZERO<RealT>;
      upper[0]              = in_service ? Pmax_ : ZERO<RealT>;
      upper[1]              = in_service ? Qmax_ : ZERO<RealT>;
      return 0;
    }

    /**
     * @brief Start from zero power
     */
    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::initialize()
    {
      ScalarT* x = x_.getData();
      x[0]       = ZERO<RealT>;
      x[1]       = ZERO<RealT>;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] typename Generator<scalar_type, index_type>::ScalarT
    Generator<scalar_type, index_type>::objective(const ScalarT* x) const
    {
      const ScalarT p = x[0];
      return in_service_ * (c0_ + c1_ * p + c2_ * p * p);
    }

    /**
     * @brief Active and reactive power into the bus
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Generator<scalar_type, index_type>::externalEqualities(const ScalarT* x, ScalarT* c) const
    {
      c[0] = x[0];
      c[1] = x[1];
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit

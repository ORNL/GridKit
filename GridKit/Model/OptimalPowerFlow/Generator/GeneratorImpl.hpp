/**
 * @file GeneratorImpl.hpp
 * @brief Definition of an optimal power flow generator.
 */

#pragma once

#include <cmath>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModelImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Generator<scalar_type, index_type>::Generator(const ModelDataT& data)
      : ComponentModelT(data.id, busNumbers(data), data.in_service),
        Pmin_(parameter(data, GeneratorParameters::Pmin, -UNBOUNDED<RealT>)),
        Pmax_(parameter(data, GeneratorParameters::Pmax, UNBOUNDED<RealT>)),
        Qmin_(parameter(data, GeneratorParameters::Qmin, -UNBOUNDED<RealT>)),
        Qmax_(parameter(data, GeneratorParameters::Qmax, UNBOUNDED<RealT>)),
        c0_(parameter(data, GeneratorParameters::c0, ZERO<RealT>)),
        c1_(parameter(data, GeneratorParameters::c1, ZERO<RealT>)),
        c2_(parameter(data, GeneratorParameters::c2, ZERO<RealT>))
    {
      // Out of service, the internal constraints fix the power at zero
      variable_lower_ = {ZERO<RealT>, ZERO<RealT>};
      variable_upper_ = {ZERO<RealT>, ZERO<RealT>};
      if (data.in_service)
      {
        variable_lower_ = {Pmin_, Qmin_};
        variable_upper_ = {Pmax_, Qmax_};
      }
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
     * @brief Start from zero power
     */
    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::initialize(ScalarT* x)
    {
      x[variable_indices_[0]] = ZERO<RealT>;
      x[variable_indices_[1]] = ZERO<RealT>;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline typename Generator<scalar_type, index_type>::ScalarT
    Generator<scalar_type, index_type>::objective(const ScalarT* x) const
    {
      const ScalarT p = x[0];
      return in_service_ * (c0_ + c1_ * p + c2_ * p * p);
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void
    Generator<scalar_type, index_type>::constraints(const ScalarT* x, ScalarT* g) const
    {
      g[0] = x[0];
      g[1] = x[1];
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit

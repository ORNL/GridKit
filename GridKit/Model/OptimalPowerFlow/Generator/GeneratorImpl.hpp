/**
 * @file GeneratorImpl.hpp
 * @brief Definition of an optimal power flow generator.
 */

#pragma once

#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Generator<scalar_type, index_type>::Generator(const ModelDataT& data)
      : ComponentT(std::type_identity<Generator>{}, data.id)
    {
      initializeParameters(data);
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Generator<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Generator");
      checks.check(Pmin_ <= Pmax_, "Pmin must not exceed Pmax");
      checks.check(Qmin_ <= Qmax_, "Qmin must not exceed Qmax");
      return checks;
    }

    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::initialize()
    {
      const auto P = static_cast<size_t>(GeneratorInternalVariables::P);
      const auto Q = static_cast<size_t>(GeneratorInternalVariables::Q);

      RealT* x_lower = x_lower_.getData();
      RealT* x_upper = x_upper_.getData();

      x_lower[P] = ZERO<RealT>;
      x_upper[P] = ZERO<RealT>;
      x_lower[Q] = ZERO<RealT>;
      x_upper[Q] = ZERO<RealT>;
      if (in_service_)
      {
        x_lower[P] = Pmin_;
        x_upper[P] = Pmax_;
        x_lower[Q] = Qmin_;
        x_upper[Q] = Qmax_;
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] typename Generator<scalar_type, index_type>::ScalarT
    Generator<scalar_type, index_type>::objective(const ScalarT* x, [[maybe_unused]] const ScalarT* external) const
    {
      const auto P = static_cast<size_t>(GeneratorInternalVariables::P);

      return u_ * (c0_ + c1_ * x[P] + c2_ * x[P] * x[P]);
    }

    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Generator<scalar_type, index_type>::externalConstraints(const ScalarT*                  x,
                                                            [[maybe_unused]] const ScalarT* external,
                                                            ScalarT*                        g) const
    {
      const auto P  = static_cast<size_t>(GeneratorInternalVariables::P);
      const auto Q  = static_cast<size_t>(GeneratorInternalVariables::Q);
      const auto DP = static_cast<size_t>(GeneratorExternalConstraints::DP);
      const auto DQ = static_cast<size_t>(GeneratorExternalConstraints::DQ);

      g[DP] = x[P];
      g[DQ] = x[Q];
    }

    template <typename scalar_type, typename index_type>
    void Generator<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = GeneratorParameters;

      const Model::ParameterReader reader(data, "Generator");

      Pmin_       = reader.readReal(Params::Pmin).value_or(Pmin_);
      Pmax_       = reader.readReal(Params::Pmax).value_or(Pmax_);
      Qmin_       = reader.readReal(Params::Qmin).value_or(Qmin_);
      Qmax_       = reader.readReal(Params::Qmax).value_or(Qmax_);
      c0_         = reader.readReal(Params::c0).value_or(c0_);
      c1_         = reader.readReal(Params::c1).value_or(c1_);
      c2_         = reader.readReal(Params::c2).value_or(c2_);
      in_service_ = reader.readBool(Params::in_service).value_or(in_service_);
      setDerivedParameters();
    }

    template <typename scalar_type, typename index_type>
    void Generator<scalar_type, index_type>::setDerivedParameters()
    {
      u_ = ZERO<RealT>;
      if (in_service_)
      {
        u_ = ONE<RealT>;
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit

/**
 * @file SignalImpl.hpp
 * @brief Definition of an optimal power flow signal.
 */

#pragma once

#include <format>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Signal<scalar_type, index_type>::Signal(IdxT signal_id)
      : ComponentT(std::type_identity<Signal>{}, std::format("signal {}", signal_id))
    {
    }

    /// A signal has no parameters
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Signal<scalar_type, index_type>::verify() const
    {
      return Model::ConfigurationChecks("Signal");
    }

    template <typename scalar_type, typename index_type>
    int Signal<scalar_type, index_type>::initialize()
    {
      const auto S = static_cast<size_t>(SignalInternalVariables::S);

      ScalarT* x = x_.getData();

      x[S] = ONE<RealT>;
      return 0;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit

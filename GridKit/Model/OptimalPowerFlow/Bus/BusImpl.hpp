/**
 * @file BusImpl.hpp
 * @brief Definition of an optimal power flow bus.
 */

#pragma once

#include <format>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : ComponentT(std::type_identity<Bus>{}, std::format("bus {}", data.number))
    {
      initializeParameters(data);
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Bus<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Bus");
      checks.check(Vmin_ >= ZERO<RealT>, "Vmin must be nonnegative");
      checks.check(Vmax_ >= Vmin_, "Vmax must not be less than Vmin");
      return checks;
    }

    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::tagInequality()
    {
      const auto P   = static_cast<size_t>(BusInternalConstraints::P);
      const auto Q   = static_cast<size_t>(BusInternalConstraints::Q);
      const auto VM2 = static_cast<size_t>(BusInternalConstraints::VM2);

      tag_[P]   = false;
      tag_[Q]   = false;
      tag_[VM2] = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize()
    {
      const auto VR  = static_cast<size_t>(BusInternalVariables::VR);
      const auto VI  = static_cast<size_t>(BusInternalVariables::VI);
      const auto VM2 = static_cast<size_t>(BusInternalConstraints::VM2);

      ScalarT* x       = x_.getData();
      RealT*   x_lower = x_lower_.getData();
      RealT*   x_upper = x_upper_.getData();
      RealT*   g_lower = g_lower_.getData();
      RealT*   g_upper = g_upper_.getData();

      x[VR] = ONE<RealT>;

      if (reference_)
      {
        x_lower[VI] = ZERO<RealT>;
        x_upper[VI] = ZERO<RealT>;
      }

      g_upper[VM2] = Vmax_ * Vmax_;
      if (Vmin_ > ZERO<RealT>)
      {
        g_lower[VM2] = Vmin_ * Vmin_;
      }
      return 0;
    }

    /// Shunt power into the bus, which starts its balance, and the squared voltage magnitude
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Bus<scalar_type, index_type>::internalConstraints(const ScalarT*                  x,
                                                      [[maybe_unused]] const ScalarT* external,
                                                      ScalarT*                        g) const
    {
      const auto VR  = static_cast<size_t>(BusInternalVariables::VR);
      const auto VI  = static_cast<size_t>(BusInternalVariables::VI);
      const auto P   = static_cast<size_t>(BusInternalConstraints::P);
      const auto Q   = static_cast<size_t>(BusInternalConstraints::Q);
      const auto VM2 = static_cast<size_t>(BusInternalConstraints::VM2);

      const ScalarT v2 = x[VR] * x[VR] + x[VI] * x[VI];

      g[P]   = -G_ * v2;
      g[Q]   = B_ * v2;
      g[VM2] = v2;
    }

    template <typename scalar_type, typename index_type>
    void Bus<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = BusParameters;

      const Model::ParameterReader reader(data, "Bus");

      Vmin_      = reader.readReal(Params::Vmin).value_or(Vmin_);
      Vmax_      = reader.readReal(Params::Vmax).value_or(Vmax_);
      G_         = reader.readReal(Params::G).value_or(G_);
      B_         = reader.readReal(Params::B).value_or(B_);
      reference_ = reader.readBool(Params::reference).value_or(reference_);
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit

/**
 * @file Signal.hpp
 * @brief Declaration of an optimal power flow signal.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Internal variables of a `Signal`
    enum class SignalInternalVariables : size_t
    {
      S, ///< \f$s\f$ Value
    };

    /// Scalar that the devices at its ports share
    template <typename scalar_type, typename index_type>
    class Signal final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

      using ComponentT::x_;

    public:
      using InternalVariablesT   = SignalInternalVariables;
      using ExternalVariablesT   = Utilities::EmptyEnum;
      using InternalConstraintsT = Utilities::EmptyEnum;
      using ExternalConstraintsT = Utilities::EmptyEnum;

      static constexpr std::array<IdxT, 0> NUM_PORT_VARIABLES{};
      static constexpr std::array<IdxT, 0> NUM_PORT_CONSTRAINTS{};

      /// A port to the signal reads its value
      static constexpr IdxT NUM_SHARED_VARIABLES = Utilities::enum_size<SignalInternalVariables>();

      explicit Signal(IdxT signal_id);

      Model::ConfigurationChecks verify() const override;
      int                        initialize() override;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit

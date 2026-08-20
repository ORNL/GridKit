#pragma once

#include <array>
#include <concepts>

#include <GridKit/Model/PhasorDynamics/Port.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit::PhasorDynamics
{
  template <typename T>
  concept PortType = requires { typename T::ScalarT; typename T::IdxT; }
                     && std::derived_from<T, Port<typename T::ScalarT, typename T::IdxT>>;

  /// Group of ports indexed by their associated signal-variable enum.
  template <PortType port_type, Utilities::SizedEnum signal_variable_type>
  class PortGroup
  {
  public:
    using PortT           = port_type;
    using SignalVariableT = signal_variable_type;
    using ScalarT         = typename PortT::ScalarT;
    using IdxT            = typename PortT::IdxT;
    PortGroup()           = default;

    static constexpr auto size() noexcept -> std::size_t
    {
      return Utilities::enum_size<SignalVariableT>();
    }

    auto operator[](SignalVariableT variable) -> PortT&
    {
      assert(Utilities::contained_within<SignalVariableT>(variable));
      return ports_[static_cast<std::size_t>(variable)];
    }

    auto operator[](SignalVariableT variable) const -> const PortT&
    {
      assert(Utilities::contained_within<SignalVariableT>(variable));
      return ports_[static_cast<std::size_t>(variable)];
    }

    template <SignalVariableT variable>
      requires(Utilities::contained_within<SignalVariableT>(variable))
    auto port() -> PortT&
    {
      return ports_[static_cast<std::size_t>(variable)];
    }

    template <SignalVariableT variable>
      requires(Utilities::contained_within<SignalVariableT>(variable))
    auto port() const -> const PortT&
    {
      return ports_[static_cast<std::size_t>(variable)];
    }

  private:
    std::array<PortT, size()> ports_{};
  };
} // namespace GridKit::PhasorDynamics

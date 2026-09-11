/**
 * @file PortGroup.hpp
 * @author superwhiskers <whiskerdev@protonmail.com>
 * @brief Group of ports indexed by a signal variable enumeration.
 */

#pragma once

#include <array>
#include <concepts>

#include <GridKit/Model/PhasorDynamics/Port.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// A type deriving from the `Port` class.
    template <typename T>
    concept PortType = requires { typename T::ScalarT; typename T::IdxT; }
                       && std::derived_from<T,
                                            Port<typename T::ScalarT,
                                                 typename T::IdxT>>;

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

      /// The size of this port group.
      static constexpr auto size() noexcept -> std::size_t
      {
        return Utilities::enum_size<SignalVariableT>();
      }

      /// Overload permitting mutable indexing into this port group by signal
      /// variable.
      auto operator[](SignalVariableT variable) -> PortT&
      {
        assert(Utilities::contained_within<SignalVariableT>(variable));
        return ports_[static_cast<std::size_t>(variable)];
      }

      /// Overload permitting indexing into this port group by signal variable.
      auto operator[](SignalVariableT variable) const -> const PortT&
      {
        assert(Utilities::contained_within<SignalVariableT>(variable));
        return ports_[static_cast<std::size_t>(variable)];
      }

      /// Checked mutable indexing into this port group by signal variable.
      template <SignalVariableT variable>
        requires(Utilities::contained_within<SignalVariableT>(variable))
      auto port() -> PortT&
      {
        return ports_[static_cast<std::size_t>(variable)];
      }

      /// Checked indexing into this port group by signal variable.
      template <SignalVariableT variable>
        requires(Utilities::contained_within<SignalVariableT>(variable))
      auto port() const -> const PortT&
      {
        return ports_[static_cast<std::size_t>(variable)];
      }

    private:
      /// The ports contained in this port group.
      std::array<PortT, size()> ports_{};
    };
  } // namespace PhasorDynamics
} // namespace GridKit

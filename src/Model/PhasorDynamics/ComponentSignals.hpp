#pragma once

#include <array>
#include <optional>
#include <type_traits>

#include <Model/PhasorDynamics/SignalNode/SignalNode.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Dummy `Variables` type for components with no variables
    enum class NoVariables : size_t
    {
      MAXIMUM
    };

    /// Concept requiring an enum to have a `MAXIMUM` variant and that it have
    /// the underlying type of `size_t`. This does not ensure the variant is
    /// the actual maximum
    template <typename T>
    concept EnumHasMaximumValueAndIsSizeT = std::is_enum_v<T>
                                            && std::is_same_v<std::underlying_type_t<T>, size_t>
                                            && requires { T::MAXIMUM; };

    /// Extension object for `Component`s adding methods and member variables
    /// related to signal bus management
    ///
    /// This is used by adding an instance in a field to your class and
    /// exposing this field to others
    template <class ScalarT, typename IdxT, typename InternalVariables, typename ExternalVariables>
      requires EnumHasMaximumValueAndIsSizeT<InternalVariables>
               && EnumHasMaximumValueAndIsSizeT<ExternalVariables>
    class ComponentSignals
    {
    public:
      /// Attaches a signal node to an external variable on this component
      template <ExternalVariables variable>
      auto attachSignalNode(SignalNode<ScalarT, IdxT>* signal)
      {
        external_variable_signals_[static_cast<size_t>(variable)] = signal;
      }

      /// Check if a signal node has been attached to an external variable
      template <ExternalVariables variable>
      auto isAttached() -> bool
      {
        return static_cast<bool>(external_variable_signals_[static_cast<size_t>(variable)]);
      }

      /// Check if a signal node has been assigned to an internal variable
      template <InternalVariables variable>
      auto isAssigned() -> bool
      {
        return static_cast<bool>(internal_variable_signals_[static_cast<size_t>(variable)]);
      }

      /// Returns a signal node for an internal signal variable to be
      /// attached to an external variable on another component
      template <InternalVariables variable>
      auto getSignalNode() -> SignalNode<ScalarT, IdxT>*
      {
        if (!internal_variable_signals_[static_cast<size_t>(variable)])
        {
          throw "A signal node has not been assigned to this internal variable";
        }

        return *internal_variable_signals_[static_cast<size_t>(variable)];
      }

      /// Returns the value of the specified external variable
      template <ExternalVariables variable>
      auto readExternalVariable() const -> ScalarT
      {
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw "A signal node has not been assigned to this external variable";
        }

        return (*external_variable_signals_[static_cast<size_t>(variable)])->read();
      }

      /// Writes a value to the specified external variable
      template <ExternalVariables variable>
      auto writeExternalVariable(ScalarT value)
      {
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw "A signal node has not been assigned to this external variable";
        }

        (*external_variable_signals_[static_cast<size_t>(variable)])->init(value);
      }

      /// Assign a signal node to an internal variable
      template <InternalVariables variable>
      auto assignSignalNode(SignalNode<ScalarT, IdxT>* node)
      {
        internal_variable_signals_[static_cast<size_t>(variable)] = node;
      }

    private:
      /// Internal variables which may have a signal associated with them for
      /// use elsewhere
      std::array<std::optional<SignalNode<ScalarT, IdxT>*>,
                 static_cast<size_t>(InternalVariables::MAXIMUM)>
          internal_variable_signals_;

      /// External variables which may have a signal associated with them for
      /// use internally
      std::array<std::optional<SignalNode<ScalarT, IdxT>*>,
                 static_cast<size_t>(ExternalVariables::MAXIMUM)>
          external_variable_signals_;
    };
  } // namespace PhasorDynamics
} // namespace GridKit

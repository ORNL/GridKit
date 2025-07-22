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

    /// "Extension" class for `Component`s adding methods and member variables
    /// related to signal bus management
    ///
    /// This was implemented in a separate class to avoid requiring users to
    /// declare signal variable enumerations when declaring a `real_type` type
    /// alias, among other things
    template <class ScalarT, typename IdxT, typename InternalVariables, typename ExternalVariables>
      requires EnumHasMaximumValueAndIsSizeT<InternalVariables>
               && EnumHasMaximumValueAndIsSizeT<ExternalVariables>
    class ComponentSignalExtension
    {
    public:
      /// Attaches a signal node to an external variable on this component
      template <ExternalVariables variable>
      auto attachSignalNode(const SignalNode<ScalarT, IdxT>& signal)
      {
        external_variable_signals_[static_cast<size_t>(variable)] = signal;
      }

      /// Returns a signal node for attaching to an external variable on
      /// another component
      template <InternalVariables variable>
      auto getSignalNode() -> const SignalNode<ScalarT, IdxT>&
      {
        if (!internal_variable_signals_[static_cast<size_t>(variable)])
        {
          throw "A signal node has not been assigned to this internal variable";
        }

        return *internal_variable_signals_[static_cast<size_t>(variable)];
      }

      /// Assigns a signal node to an internal variable
      template <InternalVariables variable>
      auto assignSignalVariable(const SignalNodeData<ScalarT, IdxT>& data)
      {
        internal_variable_signals_[static_cast<size_t>(variable)] = SignalNode(data);
      }

    private:
      /// Internal variables which may have a signal associated with them for
      /// use elsewhere
      std::array<std::optional<SignalNode<ScalarT, IdxT>>,
                 static_cast<size_t>(InternalVariables::MAXIMUM)>
          internal_variable_signals_;

      /// External variables which may have a signal associated with them for
      /// use internally
      std::array<std::optional<const SignalNode<ScalarT, IdxT>&>,
                 static_cast<size_t>(ExternalVariables::MAXIMUM)>
          external_variable_signals_;
    };
  } // namespace PhasorDynamics
} // namespace GridKit

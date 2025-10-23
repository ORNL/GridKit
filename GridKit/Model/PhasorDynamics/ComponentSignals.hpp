#pragma once

#include <array>
#include <optional>
#include <stdexcept>
#include <type_traits>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>

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
    ///
    /// @tparam ScalarT Scalar value type
    /// @tparam IdxT Index type
    /// @tparam InternalVariables An enumeration satisfying
    ///         `EnumHasMaximumValueAndIsSizeT` enumerating internal variables
    ///         for the component
    /// @tparam ExternalVariables An enumeration satisfying
    ///         `EnumHasMaximumValueAndIsSizeT` enumerating external variables
    ///         for the component
    /// @invariant InternalVariables::MAXIMUM is the greatest attainable
    ///            integer value of the enum
    /// @invariant ExternalVariables::MAXIMUM is the greatest attainable
    ///            integer value of the enum
    template <class ScalarT, typename IdxT, typename InternalVariables, typename ExternalVariables>
      requires EnumHasMaximumValueAndIsSizeT<InternalVariables>
               && EnumHasMaximumValueAndIsSizeT<ExternalVariables>
    class ComponentSignals
    {
    public:
      /// Attaches a signal node to an external variable on this component
      ///
      /// @tparam variable The external variable to attach the provided
      ///         signal to
      /// @param[in] node The signal node to attach
      /// @pre The provided pointer to a signal node is not `nullptr`
      /// @post The provided signal node is attached to the indicated
      ///       external variable
      template <ExternalVariables variable>
      auto attachSignalNode(SignalNode<ScalarT, IdxT>* node)
      {
#ifndef NDEBUG
        if (node == nullptr)
        {
          throw std::logic_error("A null pointer to a signal node has been passed to attachSignalNode");
        }
#endif

        static_assert(variable < ExternalVariables::MAXIMUM);
        external_variable_signals_[static_cast<size_t>(variable)] = node;
      }

      /// Check if a signal node has been attached to an external variable
      ///
      /// @tparam variable The external variable to check
      template <ExternalVariables variable>
      auto isAttached() const -> bool
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        return static_cast<bool>(external_variable_signals_[static_cast<size_t>(variable)]);
      }

      /// Check if a signal node has been assigned to an internal variable
      ///
      /// @tparam variable The internal variable to check
      template <InternalVariables variable>
      auto isAssigned() const -> bool
      {
        static_assert(variable < InternalVariables::MAXIMUM);
        return static_cast<bool>(internal_variable_signals_[static_cast<size_t>(variable)]);
      }

      /// Returns a signal node for an internal signal variable to be
      /// attached to an external variable on another component
      ///
      /// @tparam variable The internal variable to get the assigned
      ///         signal node of
      /// @pre A signal node has been assigned to the requested internal
      ///      variable
      template <InternalVariables variable>
      auto getSignalNode() -> SignalNode<ScalarT, IdxT>*
      {
        static_assert(variable < InternalVariables::MAXIMUM);
        if (!internal_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this internal variable");
        }

        return *internal_variable_signals_[static_cast<size_t>(variable)];
      }

      /// Returns the value of the specified external variable
      ///
      /// @tparam variable The external variable to read from
      /// @pre A signal node has been assigned to the requested external
      ///      variable
      template <ExternalVariables variable>
      auto readExternalVariable() const -> ScalarT
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this external variable");
        }

        return (*external_variable_signals_[static_cast<size_t>(variable)])->read();
      }

      /// Writes a value to the specified external variable
      ///
      /// @warning This method should be used only in component initialization
      /// methods. Use only if you know what you are doing.
      ///
      /// @tparam variable The external variable to write to
      /// @param[in] value The value to write to the signal node
      /// @pre A signal node has been assigned to the requested external
      ///      variable
      /// @post The signal node of the corresponding external variable has
      ///       the given value written to it
      template <ExternalVariables variable>
      auto writeExternalVariable(ScalarT value)
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this external variable");
        }

        (*external_variable_signals_[static_cast<size_t>(variable)])->init(value);
      }

      /// Assigns a signal node to an internal variable on this component
      ///
      /// @tparam variable The internal variable to assign the signal node to
      /// @param[in] node The signal node to assign
      /// @pre The provided pointer to a signal node is not `nullptr`
      /// @post The provided signal node is assigned to the indicated
      ///       internal variable
      template <InternalVariables variable>
      auto assignSignalNode(SignalNode<ScalarT, IdxT>* node)
      {
#ifndef NDEBUG
        if (node == nullptr)
        {
          throw std::logic_error("A null pointer to a signal node has been passed to assignSignalNode");
        }
#endif

        static_assert(variable < InternalVariables::MAXIMUM);
        internal_variable_signals_[static_cast<size_t>(variable)] = node;
      }

    private:
      /// Internal variables which may have a signal associated with them for
      /// use elsewhere
      std::array<std::optional<SignalNode<ScalarT, IdxT>*>,
                 static_cast<size_t>(InternalVariables::MAXIMUM)>
          internal_variable_signals_{};

      /// External variables which may have a signal associated with them for
      /// use internally
      std::array<std::optional<SignalNode<ScalarT, IdxT>*>,
                 static_cast<size_t>(ExternalVariables::MAXIMUM)>
          external_variable_signals_{};
    };
  } // namespace PhasorDynamics
} // namespace GridKit

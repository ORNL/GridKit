#pragma once

#include <array>
#include <optional>
#include <stdexcept>
#include <type_traits>

#include <GridKit/Model/EMT/SignalNode/SignalNode.hpp>

namespace GridKit
{
  namespace EMT
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
    /// @tparam scalar_type Scalar value type
    /// @tparam index_type Index type
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
    template <typename scalar_type, typename index_type, typename InternalVariables, typename ExternalVariables>
      requires EnumHasMaximumValueAndIsSizeT<InternalVariables>
               && EnumHasMaximumValueAndIsSizeT<ExternalVariables>
    class ComponentSignals
    {
    public:
      /// Scalar value type
      using ScalarT = scalar_type;
      /// Index type
      using IdxT    = index_type;
      /// Signal node type
      using NodeT   = SignalNode<ScalarT, IdxT>;
      /// Three-phase electrical port type
      using Port3T  = Port3<ScalarT, IdxT>;

      /// Attaches a signal node to an external variable on this component
      ///
      /// @tparam variable The external variable to attach the provided
      ///         signal to
      /// @param[in] node The signal node to attach
      /// @pre The provided pointer to a signal node is not `nullptr`
      /// @post The provided signal node is attached to the indicated
      ///       external variable
      template <ExternalVariables variable>
      auto attachSignalNode(NodeT* node)
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

      /// Attaches a three-phase port to three consecutive external variables
      /// starting at the indicated variable
      ///
      /// @tparam first The external variable of the first phase
      /// @param[in] port The three-phase port to attach
      /// @pre The provided pointer to a port is not `nullptr`
      /// @post The port's phase nodes are attached to the three consecutive
      ///       external variables starting at the indicated variable
      template <ExternalVariables first>
      auto attachPort(Port3T* port)
      {
#ifndef NDEBUG
        if (port == nullptr)
        {
          throw std::logic_error("A null pointer to a port has been passed to attachPort");
        }
#endif

        static_assert(static_cast<size_t>(first) + 2 < static_cast<size_t>(ExternalVariables::MAXIMUM));
        external_variable_signals_[static_cast<size_t>(first)]     = port->a();
        external_variable_signals_[static_cast<size_t>(first) + 1] = port->b();
        external_variable_signals_[static_cast<size_t>(first) + 2] = port->c();
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

      /// Check if a signal node has been "set"
      ///
      /// @tparam variable The external variable to check
      template <ExternalVariables variable>
      auto isLinked() const -> bool
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        return external_variable_signals_[static_cast<size_t>(variable)].value()->linked();
      }

      /// Returns a signal node for an internal signal variable to be
      /// attached to an external variable on another component
      ///
      /// @tparam variable The internal variable to get the assigned
      ///         signal node of
      /// @pre A signal node has been assigned to the requested internal
      ///      variable
      template <InternalVariables variable>
      auto getSignalNode() -> NodeT*
      {
        static_assert(variable < InternalVariables::MAXIMUM);
        if (!internal_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this internal variable");
        }

        return *internal_variable_signals_[static_cast<size_t>(variable)];
      }

      /// Returns the attached signal node of the specified external variable
      ///
      /// @tparam variable The external variable to get the attached signal
      ///         node of
      /// @pre A signal node has been attached to the requested external
      ///      variable
      template <ExternalVariables variable>
      auto getAttachedSignalNode() -> NodeT*
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been attached to this external variable");
        }

        return *external_variable_signals_[static_cast<size_t>(variable)];
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

      /// Returns the derivative of the specified external variable
      ///
      /// @tparam variable The external variable to read from
      /// @pre A signal node has been assigned to the requested external
      ///      variable and the node exposes the owning variable's derivative
      template <ExternalVariables variable>
      auto readExternalVariableDerivative() const -> ScalarT
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this external variable");
        }

        return (*external_variable_signals_[static_cast<size_t>(variable)])->readDerivative();
      }

      /// Returns the global index of the specified external variable
      ///
      /// @tparam variable The external variable to read from
      /// @pre A signal node has been assigned to the requested external
      ///      variable
      template <ExternalVariables variable>
      auto readExternalVariableIndex() const -> IdxT
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this external variable");
        }

        return (*external_variable_signals_[static_cast<size_t>(variable)])->getVariableIndex();
      }

      /// Marks derivative coupling on the specified external variable
      ///
      /// A component calls this during allocation when its equations read the
      /// derivative of the external variable, so the owning component can
      /// classify the variable as differential.
      ///
      /// @tparam variable The external variable to mark
      /// @pre A signal node has been assigned to the requested external
      ///      variable
      template <ExternalVariables variable>
      auto markDerivativeCoupling()
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this external variable");
        }

        (*external_variable_signals_[static_cast<size_t>(variable)])->markDerivativeCoupling();
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
      auto assignSignalNode(NodeT* node)
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

      /// Registers every attached external variable node with a component
      ///
      /// Called from a model's allocate() so the component base class can
      /// gather external variables, derivatives, and index maps through the
      /// attached nodes.
      ///
      /// @param[in,out] component The component to register the nodes with
      template <typename ComponentT>
      auto registerExternalVariableNodes(ComponentT& component) const
      {
        for (size_t slot = 0; slot < external_variable_signals_.size(); ++slot)
        {
          if (external_variable_signals_[slot])
          {
            component.setExternalVariableNode(static_cast<IdxT>(slot), *external_variable_signals_[slot]);
          }
        }
      }

    private:
      /// Internal variables which may have a signal associated with them for
      /// use elsewhere
      std::array<std::optional<NodeT*>,
                 static_cast<size_t>(InternalVariables::MAXIMUM)>
          internal_variable_signals_{};

      /// External variables which may have a signal associated with them for
      /// use internally
      std::array<std::optional<NodeT*>,
                 static_cast<size_t>(ExternalVariables::MAXIMUM)>
          external_variable_signals_{};
    };
  } // namespace EMT
} // namespace GridKit

#pragma once

#include <array>
#include <GridKit/Utilities/Enum.hpp>
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
    };

    /// Connections carrying network values and their residual row indices.
    /// Controller inputs and outputs use SignalPorts.
    template <typename scalar_type, typename index_type, typename InternalVariables, typename ExternalVariables>
    class ComponentSignals
    {
    public:
      /// Scalar value type
      using ScalarT = scalar_type;
      /// Index type
      using IdxT    = index_type;

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

        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<ExternalVariables>());
        external_variable_signals_[static_cast<size_t>(variable)] = node;
      }

      /// Check if a signal node has been attached to an external variable
      ///
      /// @tparam variable The external variable to check
      template <ExternalVariables variable>
      auto isAttached() const -> bool
      {
        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<ExternalVariables>());
        return static_cast<bool>(external_variable_signals_[static_cast<size_t>(variable)]);
      }

      /// Check if a signal node has been assigned to an internal variable
      ///
      /// @tparam variable The internal variable to check
      template <InternalVariables variable>
      auto isAssigned() const -> bool
      {
        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<InternalVariables>());
        return static_cast<bool>(internal_variable_signals_[static_cast<size_t>(variable)]);
      }

      /// Check if a signal node has been "set"
      ///
      /// @tparam variable The external variable to check
      template <ExternalVariables variable>
      auto isLinked() const -> bool
      {
        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<ExternalVariables>());
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
      auto getSignalNode() -> SignalNode<ScalarT, IdxT>*
      {
        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<InternalVariables>());
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
        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<ExternalVariables>());
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this external variable");
        }

        return (*external_variable_signals_[static_cast<size_t>(variable)])->read();
      }

      /// Returns the global index of the specified external variable
      ///
      /// @tparam variable The external variable to read from
      /// @pre A signal node has been assigned to the requested external
      ///      variable
      template <ExternalVariables variable>
      auto readExternalVariableIndex() const -> IdxT
      {
        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<ExternalVariables>());
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this external variable");
        }

        return (*external_variable_signals_[static_cast<size_t>(variable)])->getVariableIndex();
      }

      /// Returns the global residual index of the specified external variable
      ///
      /// @tparam variable The external variable to read from
      /// @pre A signal node has been assigned to the requested external
      ///      variable
      template <ExternalVariables variable>
      auto readExternalResidualIndex() const -> IdxT
      {
        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<ExternalVariables>());
        if (!external_variable_signals_[static_cast<size_t>(variable)])
        {
          throw std::logic_error("A signal node has not been assigned to this external variable");
        }

        return (*external_variable_signals_[static_cast<size_t>(variable)])->getResidualIndex();
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
        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<ExternalVariables>());
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

        static_assert(static_cast<size_t>(variable) < Utilities::enum_size<InternalVariables>());
        internal_variable_signals_[static_cast<size_t>(variable)] = node;
      }

    private:
      /// Internal variables which may have a signal associated with them for
      /// use elsewhere
      std::array<std::optional<SignalNode<ScalarT, IdxT>*>,
                 Utilities::enum_size<InternalVariables>()>
          internal_variable_signals_{};

      /// External variables which may have a signal associated with them for
      /// use internally
      std::array<std::optional<SignalNode<ScalarT, IdxT>*>,
                 Utilities::enum_size<ExternalVariables>()>
          external_variable_signals_{};
    };
  } // namespace PhasorDynamics
} // namespace GridKit

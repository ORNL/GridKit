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
    /// Extension object for `Component`s adding methods and member variables
    /// related to signal bus management
    ///
    /// This is used by adding an instance in a field to your class and
    /// exposing this field to others
    ///
    /// @tparam scalar_type Scalar value type
    /// @tparam index_type Index type
    /// @tparam InputPorts An enumeration of signal input ports
    /// @tparam OutputPorts An enumeration of signal output ports
    template <typename scalar_type, typename index_type, typename InputPorts, typename OutputPorts>
      requires std::is_enum_v<InputPorts>
               && std::is_enum_v<OutputPorts>
    class ComponentSignals
    {
    public:
      /// Scalar value type
      using ScalarT = scalar_type;
      /// Index type
      using IdxT    = index_type;

      /// Attaches a signal node to an input port on this component
      ///
      /// @tparam port The input port to attach the provided
      ///         signal to
      /// @param[in] node The signal node to attach
      /// @pre The provided pointer to a signal node is not `nullptr`
      /// @post The provided signal node is attached to the indicated
      ///       input port
      template <InputPorts port>
      auto attachSignalNode(SignalNode<ScalarT, IdxT>* node)
      {
#ifndef NDEBUG
        if (node == nullptr)
        {
          throw std::logic_error("A null pointer to a signal node has been passed to attachSignalNode");
        }
#endif

        static_assert(port < InputPorts::SIZE);
        input_port_signals_[static_cast<size_t>(port)] = node;
      }

      /// Check if a signal node has been attached to an input port
      ///
      /// @tparam port The input port to check
      template <InputPorts port>
      auto isAttached() const -> bool
      {
        static_assert(port < InputPorts::SIZE);
        return static_cast<bool>(input_port_signals_[static_cast<size_t>(port)]);
      }

      /// Check if a signal node has been assigned to an output port
      ///
      /// @tparam port The output port to check
      template <OutputPorts port>
      auto isAssigned() const -> bool
      {
        static_assert(port < OutputPorts::SIZE);
        return static_cast<bool>(output_port_signals_[static_cast<size_t>(port)]);
      }

      /// Check if a signal node has been "set"
      ///
      /// @tparam port The input port to check
      template <InputPorts port>
      auto isLinked() const -> bool
      {
        static_assert(port < InputPorts::SIZE);
        return input_port_signals_[static_cast<size_t>(port)].value()->linked();
      }

      /// Returns a signal node for an output port to be attached to an input
      /// port on another component
      ///
      /// @tparam port The output port to get the assigned
      ///         signal node of
      /// @pre A signal node has been assigned to the requested output port
      template <OutputPorts port>
      auto getSignalNode() -> SignalNode<ScalarT, IdxT>*
      {
        static_assert(port < OutputPorts::SIZE);
        if (!output_port_signals_[static_cast<size_t>(port)])
        {
          throw std::logic_error("A signal node has not been assigned to this output port");
        }

        return *output_port_signals_[static_cast<size_t>(port)];
      }

      /// Returns the value of the specified input port
      ///
      /// @tparam port The input port to read from
      /// @pre A signal node has been attached to the requested input port
      template <InputPorts port>
      auto readExternalVariable() const -> ScalarT
      {
        static_assert(port < InputPorts::SIZE);
        if (!input_port_signals_[static_cast<size_t>(port)])
        {
          throw std::logic_error("A signal node has not been attached to this input port");
        }

        return (*input_port_signals_[static_cast<size_t>(port)])->read();
      }

      /// Returns the global index of the specified input port
      ///
      /// @tparam port The input port to read from
      /// @pre A signal node has been attached to the requested input port
      template <InputPorts port>
      auto readExternalVariableIndex() const -> IdxT
      {
        static_assert(port < InputPorts::SIZE);
        if (!input_port_signals_[static_cast<size_t>(port)])
        {
          throw std::logic_error("A signal node has not been attached to this input port");
        }

        return (*input_port_signals_[static_cast<size_t>(port)])->getVariableIndex();
      }

      /// Writes a value to the specified input port
      ///
      /// @warning This method should be used only in component initialization
      /// methods. Use only if you know what you are doing.
      ///
      /// @tparam port The input port to write to
      /// @param[in] value The value to write to the signal node
      /// @pre A signal node has been attached to the requested input port
      /// @post The signal node of the corresponding input port has
      ///       the given value written to it
      template <InputPorts port>
      auto writeExternalVariable(ScalarT value)
      {
        static_assert(port < InputPorts::SIZE);
        if (!input_port_signals_[static_cast<size_t>(port)])
        {
          throw std::logic_error("A signal node has not been attached to this input port");
        }

        (*input_port_signals_[static_cast<size_t>(port)])->init(value);
      }

      /// Assigns a signal node to an output port on this component
      ///
      /// @tparam port The output port to assign the signal node to
      /// @param[in] node The signal node to assign
      /// @pre The provided pointer to a signal node is not `nullptr`
      /// @post The provided signal node is assigned to the indicated
      ///       output port
      template <OutputPorts port>
      auto assignSignalNode(SignalNode<ScalarT, IdxT>* node)
      {
#ifndef NDEBUG
        if (node == nullptr)
        {
          throw std::logic_error("A null pointer to a signal node has been passed to assignSignalNode");
        }
#endif

        static_assert(port < OutputPorts::SIZE);
        output_port_signals_[static_cast<size_t>(port)] = node;
      }

    private:
      /// Input ports which may have a signal attached from elsewhere
      std::array<std::optional<SignalNode<ScalarT, IdxT>*>,
                 static_cast<size_t>(InputPorts::SIZE)>
          input_port_signals_{};

      /// Output ports which may have a signal associated with them for use
      /// elsewhere
      std::array<std::optional<SignalNode<ScalarT, IdxT>*>,
                 static_cast<size_t>(OutputPorts::SIZE)>
          output_port_signals_{};
    };
  } // namespace PhasorDynamics
} // namespace GridKit

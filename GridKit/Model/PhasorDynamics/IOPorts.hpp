#pragma once

#include <array>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Base class representing a port (connection to a SignalNode)
     *
     * An input or output port may be connected to one SignalNode
     */
    template <typename scalar_type, typename index_type>
    class Port
    {
    public:
      /// Scalar type used for signals
      using ScalarT     = scalar_type;
      /// Index type
      using IdxT        = index_type;
      /// Instantiation of signal node type
      using SignalNodeT = SignalNode<ScalarT, IdxT>;

      /**
       * @brief Connect a signal node to the port
       *
       * @param node SignalNode to be connected
       */
      void connect(SignalNodeT* node) noexcept
      {
        assert(node != nullptr);
        signal_node_ = node;
      }

      /**
       * @brief Check if the port has a signal node connected
       */
      operator bool() const noexcept
      {
        return connected();
      }

      /**
       * @brief Check if the port has a signal node connected
       */
      bool connected() const noexcept
      {
        return (signal_node_ != nullptr);
      }

      /**
       * @brief Check that the connected signal node is linked with a source
       * signal
       *
       * @pre Port must be connected() to a signal node
       */
      bool linked() const noexcept
      {
        assert(connected());
        return signal_node_->linked();
      }

    protected:
      virtual void assign([[maybe_unused]] SignalNodeT*) const
      {
      }

      SignalNodeT* signal_node_{nullptr};
    };

    /**
     * @brief Port for sending a signal, making it available via a SignalNode
     */
    template <typename scalar_type, typename index_type>
    class OutputPort : public Port<scalar_type, index_type>
    {
    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using SignalNodeT = SignalNode<ScalarT, IdxT>;

      /**
       * @brief Link the connected signal node with the actual signal
       *
       * @pre Port must be connected() to a signal node
       *
       * @param signal_var Pointer to variable to access as signal
       * @param global_index Pointer to index variable providing the global
       * system index for the signal variable
       */
      void link(ScalarT* signal_var, IdxT* global_index)
      {
        assert(this->connected());
        this->signal_node_->link(signal_var, global_index);
      }

    protected:
      void assign(SignalNodeT* node) const override
      {
        node->setAssigned();
      }
    };

    /**
     * @brief Port for receiving a signal from a SignalNode
     */
    template <typename scalar_type, typename index_type>
    class InputPort : public Port<scalar_type, index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;

      /**
       * @brief Read the signal (value) from the connected SignalNode
       *
       * @pre Port must be connected() to a signal node
       */
      ScalarT readSignal() const
      {
        assert(this->connected());
        return this->signal_node_->read();
      }

      /**
       * @brief Get the global index of the variable provided as signal
       */
      IdxT signalVariableIndex() const
      {
        assert(*this);
        return this->signal_node_->getVariableIndex();
      }

      /**
       * @brief Write a value back to the sending port
       *
       * @warning This is hazardous and should at most be used during
       * initialization methods.
       *
       * @todo This is currently required by Genrou and Gensal components (or
       * perhaps more accurately, for Tgov1 to be properly initialized), but
       * it thoroughly violates the concept being implemented. We need to find a
       * different way for those component types to handle the synchronization
       * required at initialization.
       */
      void writeValue(ScalarT value)
      {
        this->signal_node_->init(value);
      }
    };

    /**
     * @brief Represents a group of related ports
     *
     * This is a grouping of input or output ports. Each port is associated with
     * a specific signal variable that is expected to be passed through the
     * signal node connected to that port.
     *
     * @tparam port_type InputPort or OutputPort (with template arguments)
     * @tparam An enumeration that specifies the group of ports (must satisfy
     * the SizedEnum concept)
     */
    template <typename port_type, Utilities::SizedEnum port_enum_type>
    class PortGroup
    {
    public:
      using PortT          = port_type;
      using PortEnumT      = port_enum_type;
      using ScalarT        = PortT::ScalarT;
      using IdxT           = PortT::IdxT;
      using SignalNodeT    = SignalNode<ScalarT, IdxT>;
      using SignalNodeSetT = SignalNodeSet<SignalNodeT>;

      PortGroup() = default;

      /**
       * @brief Connect signal nodes mapped from signal ids provided in map
       *
       * @param id_map Map of PortEnumT variables to signal ids
       * @param signal_nodes SignalNodeSet mapping signal ids to signal nodes
       */
      template <typename signal_id_map_type>
      PortGroup(const signal_id_map_type& id_map, SignalNodeSetT& signal_nodes)
      {
        for (std::size_t i = 0; i < Utilities::enum_size<PortEnumT>(); ++i)
        {
          auto var = static_cast<PortEnumT>(i);
          if (id_map.contains(var))
          {
            (*this)[var].connect(signal_nodes[id_map.at(var)]);
          }
        }
      }

      /**
       * @brief Get number or ports
       */
      static constexpr std::size_t size() noexcept
      {
        return Utilities::enum_size<PortEnumT>();
      }

      ///@{
      /**
       * @brief Access port vor specified signal variable
       */
      PortT& operator[](PortEnumT var)
      {
        return ports_[static_cast<std::size_t>(var)];
      }

      const PortT& operator[](PortEnumT var) const
      {
        return ports_[static_cast<std::size_t>(var)];
      }

      ///@}

    private:
      std::array<PortT, size()> ports_{};
    };

    /**
     * @brief Set of input and output ports for a component
     *
     * @tparam scalar_type Scalar value type used for signals
     * @tparam model_data_type "<model>Data" type specifying variables. This is
     * expected to be (or derive from) an instantiation of ComponentData. This
     * provides the SignalInputs and SignalOutputs enumeration types for
     * specifying each PortGroup
     */
    template <typename scalar_type, typename model_data_type>
    struct IOPorts
    {
      /// Scalar type used for signals
      using ScalarT        = scalar_type;
      /// Model data specification type
      using ModelDataT     = model_data_type;
      /// Index type
      using IdxT           = ModelDataT::IdxT;
      /// Enum specifying signal variables for input
      using SignalInEnumT  = ModelDataT::SignalInputs;
      /// Enum specifying signal variables for output
      using SignalOutEnumT = ModelDataT::SignalOutputs;
      /// Input port type
      using InputPortT     = InputPort<ScalarT, IdxT>;
      /// Output port type
      using OutputPortT    = OutputPort<ScalarT, IdxT>;
      /// SignalNode type
      using SignalNodeT    = SignalNode<ScalarT, IdxT>;
      /// SignalNodeSet type
      using SignalNodeSetT = SignalNodeSet<SignalNodeT>;

      /// Grouping of input ports
      PortGroup<InputPortT, SignalInEnumT>   in;
      /// Grouping of output ports
      PortGroup<OutputPortT, SignalOutEnumT> out;

      IOPorts() = default;

      /**
       * @brief Construct with model data and SignalNodeSet
       *
       * This will connect signal nodes to input and output ports according to
       * the signal ids specified in the model data's `signal_inputs` and
       * `signal_outputs` id maps.
       */
      IOPorts(const ModelDataT& data, SignalNodeSetT& signal_nodes)
        : in{data.signal_inputs, signal_nodes},
          out{data.signal_outputs, signal_nodes}
      {
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit

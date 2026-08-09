/**
 * @file SignalNode model implementation.
 */
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    SignalNode<scalar_type, index_type>::SignalNode()
    {
    }

    template <typename scalar_type, typename index_type>
    SignalNode<scalar_type, index_type>::SignalNode(const SignalNodeData<RealT, IdxT>& data)
      : signal_id_(data.signal_id)
    {
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::set(ScalarT* signal, IdxT* variable_index)
    {
      signal_         = signal;
      variable_index_ = variable_index;
    }

    template <typename scalar_type, typename index_type>
    bool SignalNode<scalar_type, index_type>::linked() const
    {
      return (signal_) && (variable_index_);
    }

    template <typename scalar_type, typename index_type>
    scalar_type SignalNode<scalar_type, index_type>::read() const
    {
      return *signal_;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::init(ScalarT signal)
    {
      if (isJunction() && initializing_)
      {
        throw std::logic_error("A cycle was encountered while initializing a signal-node junction");
      }

      *signal_ = signal;

      if (!isJunction())
      {
        return;
      }

      initializing_ = true;
      try
      {
        const auto initialization_input = static_cast<std::size_t>(junction_initialization_input_index_);

        ScalarT remainder{junction_bias_};
        for (std::size_t i = 0; i < junction_inputs_.size(); ++i)
        {
          if (i != initialization_input)
          {
            const auto& input  = junction_inputs_[i];
            remainder         += input.gain * input.signal->read();
          }
        }

        const auto& input = junction_inputs_[initialization_input];
        input.signal->init((signal - remainder) / input.gain);
      }
      catch (...)
      {
        initializing_ = false;
        throw;
      }
      initializing_ = false;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::configureJunction(
        RealT                      bias,
        IdxT                       initialization_input_index,
        std::vector<JunctionInput> inputs)
    {
      if (isJunction())
      {
        throw std::logic_error("A signal-node junction cannot be configured more than once");
      }
      if (inputs.empty())
      {
        throw std::invalid_argument("A signal-node junction requires at least one input");
      }
      if (!std::isfinite(bias))
      {
        throw std::invalid_argument("Signal-node junction bias must be finite");
      }
      if constexpr (std::is_signed_v<IdxT>)
      {
        if (initialization_input_index < 0)
        {
          throw std::invalid_argument("Signal-node junction initialization input is out of range");
        }
      }

      const auto initialization_input = static_cast<std::size_t>(initialization_input_index);
      if (initialization_input >= inputs.size())
      {
        throw std::invalid_argument("Signal-node junction initialization input is out of range");
      }

      for (std::size_t i = 0; i < inputs.size(); ++i)
      {
        if (inputs[i].signal == nullptr)
        {
          throw std::invalid_argument("A signal-node junction input cannot be null");
        }
        if (inputs[i].signal == this)
        {
          throw std::invalid_argument("A signal-node junction cannot directly reference itself");
        }
        if (!std::isfinite(inputs[i].gain))
        {
          throw std::invalid_argument("Signal-node junction input gains must be finite");
        }
        for (std::size_t j = 0; j < i; ++j)
        {
          if (inputs[i].signal == inputs[j].signal)
          {
            throw std::invalid_argument("A signal node cannot occur more than once in a junction");
          }
        }
      }

      if (inputs[initialization_input].gain == ZERO<RealT>)
      {
        throw std::invalid_argument("Signal-node junction initialization input gain cannot be zero");
      }

      junction_bias_                       = bias;
      junction_initialization_input_index_ = initialization_input_index;
      junction_inputs_                     = std::move(inputs);
    }

    template <typename scalar_type, typename index_type>
    bool SignalNode<scalar_type, index_type>::isJunction() const
    {
      return !junction_inputs_.empty();
    }

    template <typename scalar_type, typename index_type>
    typename SignalNode<scalar_type, index_type>::RealT
    SignalNode<scalar_type, index_type>::junctionBias() const
    {
      return junction_bias_;
    }

    template <typename scalar_type, typename index_type>
    typename SignalNode<scalar_type, index_type>::IdxT
    SignalNode<scalar_type, index_type>::junctionInitializationInputIndex() const
    {
      return junction_initialization_input_index_;
    }

    template <typename scalar_type, typename index_type>
    const std::vector<typename SignalNode<scalar_type, index_type>::JunctionInput>&
    SignalNode<scalar_type, index_type>::junctionInputs() const
    {
      return junction_inputs_;
    }

    template <typename scalar_type, typename index_type>
    scalar_type SignalNode<scalar_type, index_type>::junctionValue() const
    {
      ScalarT value{junction_bias_};
      for (const auto& input : junction_inputs_)
      {
        value += input.gain * input.signal->read();
      }
      return value;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::initializeJunction()
    {
      if (!isJunction())
      {
        throw std::logic_error("Cannot initialize an ordinary signal node as a junction");
      }
      *signal_ = junctionValue();
    }
  } // namespace PhasorDynamics
} // namespace GridKit

#pragma once

#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename real_type, typename index_type>
    struct SignalNodeData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief SignalNode model implementation base class.
     *
     */
    template <typename scalar_type, typename index_type>
    class SignalNode
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      /// One resolved input to an algebraic signal-node junction.
      struct JunctionInput
      {
        SignalNode* signal{nullptr};
        RealT       gain{ONE<RealT>};
      };

      SignalNode();
      SignalNode(const SignalNodeData<RealT, IdxT>& data);

      virtual ~SignalNode() = default;

      void    set(ScalarT* signal_in, IdxT* global_index);
      bool    linked() const;
      ScalarT read() const;
      void    init(ScalarT signal_in);

      /**
       * @brief Configure this signal node as an algebraic junction.
       *
       * @param[in] bias Constant offset in the junction equation.
       * @param[in] initialization_input_index Resolved position of the input
       *            to back-solve when @ref init writes a requested output.
       * @param[in] inputs Resolved, weighted input signal nodes.
       *
       * @pre The input list is nonempty, contains no null, duplicate, or
       *      direct self references, and the initialization input has a
       *      nonzero gain.
       */
      void configureJunction(RealT                      bias,
                             IdxT                       initialization_input_index,
                             std::vector<JunctionInput> inputs);

      bool                              isJunction() const;
      RealT                             junctionBias() const;
      IdxT                              junctionInitializationInputIndex() const;
      const std::vector<JunctionInput>& junctionInputs() const;
      ScalarT                           junctionValue() const;
      void                              initializeJunction();

      const IdxT signalId() const
      {
        return signal_id_;
      }

      IdxT getVariableIndex() const
      {
        return *variable_index_;
      }

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

    private:
      ScalarT* signal_{nullptr};
      IdxT     signal_id_{0};

      RealT                      junction_bias_{ZERO<RealT>};
      IdxT                       junction_initialization_input_index_{0};
      std::vector<JunctionInput> junction_inputs_;
      bool                       initializing_{false};

    protected:
      const IdxT bus_id_{INVALID_INDEX<IdxT>};

      IdxT* variable_index_{nullptr};
    };

  } // namespace PhasorDynamics
} // namespace GridKit

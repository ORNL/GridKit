/**
 * @file SignalNodeData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for signal nodes
 *
 */
#pragma once

#include <optional>
#include <string>
#include <vector>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// One weighted input to a signal-node junction.
    template <typename real_type, typename index_type>
    struct SignalNodeJunctionInputData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      IdxT  signal_id{0};   ///< Identifier of the input signal node
      RealT gain{RealT{1}}; ///< Multiplier applied to the input signal
    };

    /**
     * @brief Modeling data for an algebraic signal-node junction.
     *
     * The junction output is constrained by
     * \f$y = b + \sum_i g_i u_i\f$. During component initialization, writes
     * to the output are propagated backward through the designated input.
     */
    template <typename real_type, typename index_type>
    struct SignalNodeJunctionData
    {
      using RealT = real_type;
      using IdxT  = index_type;
      using Input = SignalNodeJunctionInputData<RealT, IdxT>;

      RealT              bias{RealT{0}};          ///< Constant junction offset
      IdxT               initialization_input{0}; ///< Signal ID initialized from the output
      std::vector<Input> inputs;                  ///< Weighted junction inputs
    };

    /**
     * @brief Contains modeling data for a signal node.
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     */
    template <typename real_type, typename index_type>
    struct SignalNodeData
    {
      using RealT    = real_type;
      using IdxT     = index_type;
      using Junction = SignalNodeJunctionData<RealT, IdxT>;

      std::string             name;         ///< A name given to this signal node
      IdxT                    signal_id{0}; ///< The unique ID of the signal node
      std::optional<Junction> junction;     ///< Optional algebraic junction definition
    };
  } // namespace PhasorDynamics
} // namespace GridKit

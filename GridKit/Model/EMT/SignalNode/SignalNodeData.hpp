/**
 * @file SignalNodeData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for signal nodes
 *
 */
#pragma once

#include <string>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Contains modeling data for a signal node
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct SignalNodeData
    {
      using IdxT = index_type;

      std::string name;         ///< A name given to this signal node
      IdxT        signal_id{0}; ///< The unique ID of the signal node
    };
  } // namespace EMT
} // namespace GridKit

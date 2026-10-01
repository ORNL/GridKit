/**
 * @file SignalNodeData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (nodes)
 *
 */
#pragma once

#include <bitset>
#include <optional>
#include <string>
#include <type_traits>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Contains modeling data for a Bus
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename real_type, typename index_type>
    struct SignalNodeData
    {
      using IdxT = index_type;

      std::string name;         ///< A name given to this bus
      IdxT        signal_id{0}; ///< The unique ID of the signal node
    };
  } // namespace PhasorDynamics
} // namespace GridKit

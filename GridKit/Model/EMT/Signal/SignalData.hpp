/**
 * @file SignalData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for signals
 *
 */
#pragma once

#include <optional>
#include <string>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Contains modeling data for a signal
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct SignalData
    {
      std::string              id;    ///< Unique signal identifier within the system
      std::optional<real_type> value; ///< Constant value, otherwise supplied by a producer
    };
  } // namespace EMT
} // namespace GridKit

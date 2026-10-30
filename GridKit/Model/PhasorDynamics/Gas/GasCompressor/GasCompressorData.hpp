/**
 * @file GasCompressorData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for gas compressors.
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /// Parameters for a gas compressor
      enum class GasCompressorParameters : size_t
      {
        ratio, ///< \f$\alpha\f$ Compression ratio [-]
      };

      /// Buses for a gas compressor
      enum class GasCompressorBuses : size_t
      {
        bus1, ///< \f$p_1\f$ Required Known suction pressure [Pa]
        bus2, ///< \f$p_2\f$ Required Known discharge pressure [Pa]
      };

      /// Signal inputs supported for a gas compressor
      enum class GasCompressorSignalInputs : size_t
      {
      };

      /// Signal outputs supported for a gas compressor
      enum class GasCompressorSignalOutputs : size_t
      {
      };

      /// Variables able to be monitored for a gas compressor
      enum class GasCompressorMonitorableVariables : size_t
      {
        q, ///< \f$q\f$ Mass flow from bus 1 to bus 2 [kg/s]
      };

      /**
       * @brief Contains modeling data for a GasCompressor
       *
       * @tparam real_type  Real parameter data type
       * @tparam index_type Integer parameter data type
       */
      template <typename real_type, typename index_type>
      using GasCompressorData =
          ComponentData<real_type,
                        index_type,
                        GasCompressorParameters,
                        GasCompressorBuses,
                        GasCompressorSignalInputs,
                        GasCompressorSignalOutputs,
                        GasCompressorMonitorableVariables>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit

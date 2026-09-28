/**
 * @file DependentNortonData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for a controlled Norton source.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Source
    {
      /// Required parameters for a Norton source.
      enum class DependentNortonParameters : size_t
      {
        G, ///< \f$G\f$ Parallel conductance on system base [p.u.]
        B, ///< \f$B\f$ Parallel susceptance on system base [p.u.]
      };

      /// Buses for a Norton source.
      enum class DependentNortonBuses : size_t
      {
        bus, ///< Terminal bus ID
      };

      /// Required signal inputs for a Norton source.
      enum class DependentNortonSignalInputs : size_t
      {
        inr, ///< \f$I_r^\mathrm{N}\f$ Required Known source-current input, real component on system base [p.u.]
        ini, ///< \f$I_i^\mathrm{N}\f$ Required Known source-current input, imaginary component on system base [p.u.]
      };

      /// Signal outputs for a Norton source.
      enum class DependentNortonSignalOutputs : size_t
      {
      };

      /// Variables available through the monitor interface.
      enum class DependentNortonMonitorableVariables : size_t
      {
        ir, ///< Terminal current, real component on system base [p.u.]
        ii, ///< Terminal current, imaginary component on system base [p.u.]
        p,  ///< Terminal active power on system base [p.u.]
        q,  ///< Terminal reactive power on system base [p.u.]
      };

      /**
       * @brief Parameters, terminal connection, current inputs, and monitors.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      using DependentNortonData =
          ComponentData<real_type,
                        index_type,
                        DependentNortonParameters,
                        DependentNortonBuses,
                        DependentNortonSignalInputs,
                        DependentNortonSignalOutputs,
                        DependentNortonMonitorableVariables>;
    } // namespace Source
  } // namespace PhasorDynamics
} // namespace GridKit

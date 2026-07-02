/**
 * @file Ieeet1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 * @brief Data structure for IEEET1 Data
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Initial parameters for IEEET1 Exciter model
      enum class Ieeet1Parameters
      {
        Tr,     ///< Time constant for voltage sensing
        Ka,     ///< Coefficient for voltage regulation
        Ta,     ///< Time constant for voltage regulation
        Ke,     ///< Coefficient for excitation system
        Te,     ///< Time constant for excitation system
        Kf,     ///< Coefficient for feedback
        Tf,     ///< Time constant for feedback
        Vrmin,  ///< LL to voltage regulation
        Vrmax,  ///< HH to voltage regulation
        E1,     ///< Saturation parameter
        E2,     ///< Saturation parameter
        Se1,    ///< Saturation parameter
        Se2,    ///< Saturation parameter
        Ispdlim ///< Speed limit flag indicator
      };

      /// Buses for a IEEET1 Exciter model
      enum class Ieeet1Buses : size_t
      {
        bus, ///< Unique ID of the terminal bus
        SIZE
      };

      /// Signal inputs for a IEEET1 Exciter model
      enum class Ieeet1SignalInputs : size_t
      {
        speed, ///< Unique ID of the generator speed signal
        vs,    ///< Unique ID of the stabilizer output signal (optional)
        SIZE
      };

      /// Signal outputs for a IEEET1 Exciter model
      enum class Ieeet1SignalOutputs : size_t
      {
        efd, ///< Unique ID of the output efd signal
        SIZE
      };

      /// Variables able to be monitored for a IEEET1 Exciter model
      enum class Ieeet1MonitorableVariables
      {
        efd,
        ksat
      };

      /**
       * @brief Contains modeling data for a IEEET1 Exciter model.
       *
       * @tparam real_type  Real parameter data type
       * @tparam index_type Integer parameter data type
       *
       * Integer parameters are of the same type as matrix and vector indices.
       */
      template <typename real_type, typename index_type>
      struct Ieeet1Data : public ComponentData<real_type,
                                               index_type,
                                               Ieeet1Parameters,
                                               Ieeet1Buses,
                                               Ieeet1SignalInputs,
                                               Ieeet1SignalOutputs,
                                               Ieeet1MonitorableVariables>
      {
        Ieeet1Data() = default;

        using Parameters           = Ieeet1Parameters;
        using Buses                = Ieeet1Buses;
        using SignalInputs         = Ieeet1SignalInputs;
        using SignalOutputs        = Ieeet1SignalOutputs;
        using MonitorableVariables = Ieeet1MonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

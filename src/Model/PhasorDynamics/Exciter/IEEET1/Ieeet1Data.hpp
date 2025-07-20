/**
 * @file Ieeet1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * 
 * @brief Data structure for IEEET1 Data
 *
 */
#pragma once

#include <Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Initial parameters for a Genrou generator model
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

      /// Ports for a Genrou generator model
      enum class Ieeet1Ports
      {
        bus,           ///< Unique ID of the terminal bus
        speed_signal,  ///< Unique ID of the generator speed signal
        efd_signal,    ///< Unique ID of the output efd signal
      };

      /// Variables able to be monitored for a Genrou generator model
      enum class Ieeet1MonitorableVariables
      {
        efd,
        ksat,
      };

      /**
       * @brief Contains modeling data for a Genrou generator model.
       *
       * @tparam RealT Real parameter data type
       * @tparam IdxT  Integer parameter data type
       *
       * Integer parameters are of the same type as matrix and vector indices.
       */
      template <typename RealT, typename IdxT>
      struct Ieeet1Data : public ComponentData<RealT,
                                              IdxT,
                                              Ieeet1Parameters,
                                              Ieeet1Ports,
                                              Ieeet1MonitorableVariables>
      {
        Ieeet1Data() = default;

        using Parameters           = Ieeet1Parameters;
        using Ports                = Ieeet1Ports;
        using MonitorableVariables = Ieeet1MonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

/**
 * @file Ieeet1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 * @brief Data structure for IEEET1 Data
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

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

      /// Ports for a IEEET1 Exciter model
      enum class Ieeet1Ports
      {
        bus,   ///< Unique ID of the terminal bus
        speed, ///< Unique ID of the generator speed signal
        efd,   ///< Unique ID of the output efd signal
      };

      /// Variables able to be monitored for a IEEET1 Exciter model
      enum class Ieeet1MonitorableVariables
      {
        efd,
        ksat,
        SIZE
      };

      /**
       * @brief Convert enum value to string label
       */
      inline const std::string& enumLabel(Ieeet1MonitorableVariables v)
      {
        static std::string labels[] = {"efd", "ksat"};
        return labels[Utilities::enumId(v)];
      }

      /**
       * @brief Contains modeling data for a IEEET1 Exciter model.
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

/**
 * @file Ieeet1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 * @brief Data structure for IEEET1 Data
 *
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Initial parameters for IEEET1 Exciter model
      enum class Ieeet1Parameters
      {
        V,      ///< Rated line-to-line RMS terminal voltage in volts
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

      /// Signal inputs for a IEEET1 Exciter model
      enum class Ieeet1Inputs : size_t
      {
        bus,   ///< Terminal voltage port
        speed, ///< Unique ID of the machine rotor-speed signal (1 p.u. synchronous)
        vref,  ///< Unique ID of the voltage reference signal (optional)
        vs,    ///< Unique ID of the stabilizer output signal (optional)
        vuel,  ///< Unique ID of the under-excitation limiter signal (optional)
        voel,  ///< Unique ID of the over-excitation limiter signal (optional)
        SIZE
      };

      /// Signal outputs for a IEEET1 Exciter model
      enum class Ieeet1Outputs : size_t
      {
        efd, ///< Unique ID of the output efd signal
        SIZE
      };

      /// Variables able to be monitored for a IEEET1 Exciter model
      enum class Ieeet1MonitorableVariables
      {
        efd,
        ksat,
        vts,
        vr,
        vref
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
                                               Ieeet1Inputs,
                                               Ieeet1Outputs,
                                               Ieeet1MonitorableVariables>
      {
        Ieeet1Data() = default;

        using Parameters           = Ieeet1Parameters;
        using Inputs               = Ieeet1Inputs;
        using Outputs              = Ieeet1Outputs;
        using MonitorableVariables = Ieeet1MonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

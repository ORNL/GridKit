/**
 * @file Ieeet1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 * @brief Data structure and manifest for IEEET1 Data
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
      enum class Ieeet1Parameters : size_t
      {
        Tr,
        Ka,
        Ta,
        Ke,
        Te,
        Kf,
        Tf,
        Vrmin,
        Vrmax,
        E1,
        E2,
        Se1,
        Se2,
        Ispdlim,
      };

      /// Buses for a IEEET1 Exciter model
      enum class Ieeet1Buses : size_t
      {
        bus,
      };

      /// Signal inputs for a IEEET1 Exciter model
      enum class Ieeet1SignalInputs : size_t
      {
        speed,
        vref,
        vs,
        vuel,
        voel,
      };

      /// Signal outputs for a IEEET1 Exciter model
      enum class Ieeet1SignalOutputs : size_t
      {
        efd,
      };

      /// Variables able to be monitored for a IEEET1 Exciter model
      enum class Ieeet1MonitorableVariables : size_t
      {
        efd,
        ksat,
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
      struct Ieeet1Data : ComponentData<real_type,
                                        index_type,
                                        Ieeet1Parameters,
                                        Ieeet1Buses,
                                        Ieeet1SignalInputs,
                                        Ieeet1SignalOutputs,
                                        Ieeet1MonitorableVariables>
      {
        /// Units, symbols and descriptions of the IEEET1 data
        static constexpr typename Ieeet1Data::Manifest MANIFEST{
            .class_name  = "Ieeet1",
            .description = "IEEE Type 1 excitation-system model",
            .parameters  = {{
                {Ieeet1Parameters::Tr, "s", "T_R", "Time constant for voltage sensing"},
                {Ieeet1Parameters::Ka, "p.u.", "K_A", "Coefficient for voltage regulation"},
                {Ieeet1Parameters::Ta, "s", "T_A", "Time constant for voltage regulation"},
                {Ieeet1Parameters::Ke, "p.u.", "K_E", "Exciter field resistance line slope margin, 0 requests automatic calculation"},
                {Ieeet1Parameters::Te, "s", "T_E", "Time constant for excitation system"},
                {Ieeet1Parameters::Kf, "p.u.", "K_F", "Coefficient for feedback"},
                {Ieeet1Parameters::Tf, "s", "T_F", "Time constant for feedback"},
                {Ieeet1Parameters::Vrmin, "p.u.", "V_R^{\\min}", "Lower limit to voltage regulation"},
                {Ieeet1Parameters::Vrmax, "p.u.", "V_R^{\\max}", "Upper limit to voltage regulation"},
                {Ieeet1Parameters::E1, "p.u.", "E_1", "Saturation parameter"},
                {Ieeet1Parameters::E2, "p.u.", "E_2", "Saturation parameter"},
                {Ieeet1Parameters::Se1, "p.u.", "S_1", "Saturation parameter"},
                {Ieeet1Parameters::Se2, "p.u.", "S_2", "Saturation parameter"},
                {Ieeet1Parameters::Ispdlim, "binary", "I_{\\mathrm{spdlim}}", "Speed limit flag indicator"},
            }},
            .buses       = {{
                {Ieeet1Buses::bus, "Terminal bus voltage"},
            }},
            .inputs      = {{
                {Ieeet1SignalInputs::speed, "Machine speed deviation"},
                {Ieeet1SignalInputs::vref, "Voltage-control reference"},
                {Ieeet1SignalInputs::vs, "Stabilizer input signal"},
                {Ieeet1SignalInputs::vuel, "Under-excitation limiter input"},
                {Ieeet1SignalInputs::voel, "Over-excitation limiter input"},
            }},
            .outputs     = {{
                {Ieeet1SignalOutputs::efd, "Field-voltage output seeded by the machine"},
            }},
            .monitors    = {{
                {Ieeet1MonitorableVariables::efd, "p.u.", "Field winding voltage"},
                {Ieeet1MonitorableVariables::ksat, "p.u.", "Scaled-quadratic saturation contribution"},
            }},
        };

        static_assert(MANIFEST.complete());
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

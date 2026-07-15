/**
 * @file Esdc1aData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the ESDC1A exciter model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ESDC1A exciter model.
      enum class Esdc1aParameters
      {
        Tr,     ///< Transducer time constant
        Ka,     ///< Voltage-regulator gain
        Ta,     ///< Voltage-regulator time constant
        Tb,     ///< Lead-lag denominator time constant
        Tc,     ///< Lead-lag numerator time constant
        Vrmax,  ///< Maximum voltage-regulator output
        Vrmin,  ///< Minimum voltage-regulator output
        Ke,     ///< Exciter field-resistance line-slope margin
        Te,     ///< Exciter field time constant
        Kf,     ///< Stabilizing feedback gain
        Tf1,    ///< Feedback lead time constant
        Spdmlt, ///< Speed multiplier flag
        E1,     ///< First saturation voltage point
        Se1,    ///< Saturation value at E1
        E2,     ///< Second saturation voltage point
        Se2,    ///< Saturation value at E2
        UEL,    ///< UEL input-location selector
        exclim  ///< Exciter feedback lower-limit flag
      };

      /// Buses for the ESDC1A exciter model.
      enum class Esdc1aBuses : size_t
      {
        bus, ///< Unique ID of the terminal bus
        SIZE
      };

      /// Signal inputs for the ESDC1A exciter model.
      enum class Esdc1aSignalInputs : size_t
      {
        speed, ///< Unique ID of the generator speed-deviation signal
        vref,  ///< Unique ID of the voltage-reference signal
        vs,    ///< Unique ID of the optional stabilizer input signal
        vuel,  ///< Unique ID of the optional UEL input signal
        SIZE
      };

      /// Signal outputs for the ESDC1A exciter model.
      enum class Esdc1aSignalOutputs : size_t
      {
        efd, ///< Unique ID of the output EFD signal
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class Esdc1aMonitorableVariables
      {
        efd, ///< Field-voltage output
        vc,  ///< Sensed compensated voltage
        vr,  ///< Voltage-regulator output
        vf,  ///< Stabilizing feedback state
        se,  ///< Saturation coefficient
        vfe  ///< Exciter feedback signal
      };

      template <typename RealT, typename IdxT>
      struct Esdc1aData : public ComponentData<RealT,
                                               IdxT,
                                               Esdc1aParameters,
                                               Esdc1aBuses,
                                               Esdc1aSignalInputs,
                                               Esdc1aSignalOutputs,
                                               Esdc1aMonitorableVariables>
      {
        Esdc1aData() = default;

        using Parameters           = Esdc1aParameters;
        using Buses                = Esdc1aBuses;
        using SignalInputs         = Esdc1aSignalInputs;
        using SignalOutputs        = Esdc1aSignalOutputs;
        using MonitorableVariables = Esdc1aMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

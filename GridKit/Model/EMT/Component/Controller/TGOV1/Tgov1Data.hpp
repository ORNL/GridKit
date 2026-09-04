/**
 * @file Tgov1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the EMT TGOV1 governor
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /**
       * @brief Parameter keys for TGOV1 Governor model.
       *
       * These enum values serve as keys for the parameters map in ComponentData.
       */
      enum class Tgov1Parameters
      {
        Trate, ///< \f$T_\mathrm{rate}\f$ Turbine-rating power base
        R,     ///< \f$R\f$ Permanent droop
        T1,    ///< \f$T_1\f$ Steam-bowl time constant
        T2,    ///< \f$T_2\f$ Turbine numerator time constant
        T3,    ///< \f$T_3\f$ Reheater time constant
        Pvmax, ///< \f$P_v^\mathrm{max}\f$ Maximum valve position
        Pvmin, ///< \f$P_v^\mathrm{min}\f$ Minimum valve position
        Dt     ///< \f$D_t\f$ Turbine damping coefficient
      };

      /**
       * @brief TGOV1 inputs.
       */
      enum class Tgov1Inputs : size_t
      {
        speed, ///< \f$\omega_r\f$ Optional machine rotor-speed signal ID
        pref,  ///< \f$P_\mathrm{ref}\f$ Optional governor-reference signal ID
        SIZE,
      };

      /**
       * @brief TGOV1 outputs.
       */
      enum class Tgov1Outputs : size_t
      {
        pmech, ///< \f$P_m\f$ Required mechanical-power output signal ID
        SIZE,
      };

      /**
       * @brief Placeholder enum for TGOV1 monitorable variables.
       */
      enum class Tgov1MonitorableVariables
      {
        NONE,
      };

      /**
       * @brief Modeling data for TGOV1 Governor using ComponentData base.
       *
       * @tparam real_type  Real number type (e.g., double)
       * @tparam index_type Index type (e.g., size_t)
       */
      template <typename real_type, typename index_type>
      struct Tgov1Data : public ComponentData<real_type,
                                              index_type,
                                              Tgov1Parameters,
                                              Tgov1Inputs,
                                              Tgov1Outputs,
                                              Tgov1MonitorableVariables>
      {
        Tgov1Data() = default;

        using Parameters           = Tgov1Parameters;
        using Inputs               = Tgov1Inputs;
        using Outputs              = Tgov1Outputs;
        using MonitorableVariables = Tgov1MonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

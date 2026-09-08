/**
 * @file PllData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the EMT phase-locked loop.
 */
#pragma once
#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class PllParameters
    {
      V,  ///< Rated line-to-line RMS voltage [V]
      f,  ///< Nominal frequency [Hz]
      Kp, ///< \f$K_P\f$ Proportional gain [rad/s]
      Ki, ///< \f$K_I\f$ Integral gain [rad/s^2]
    };
    enum class PllInputs : size_t
    {
      va,
      vb,
      vc,
      SIZE
    };
    enum class PllOutputs : size_t
    {
      theta, ///< \f$\theta\f$ Electrical reference angle [rad]
      omega, ///< \f$\omega\f$ Electrical angular frequency [rad/s]
      SIZE,
    };
    enum class PllMonitorableVariables
    {
      theta,
      xi,
      omega,
      vq
    };

    template <typename real_type, typename index_type>
    struct PllData : public ComponentData<real_type, index_type, PllParameters, PllInputs, PllOutputs, PllMonitorableVariables>
    {
      PllData() = default;

      using Parameters           = PllParameters;
      using Inputs               = PllInputs;
      using Outputs              = PllOutputs;
      using MonitorableVariables = PllMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit

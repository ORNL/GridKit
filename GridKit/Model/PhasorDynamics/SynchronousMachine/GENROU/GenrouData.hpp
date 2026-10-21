/**
 * @file GenrouData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data and manifest for the GENROU synchronous machine
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a Genrou generator model
    enum class GenrouParameters : size_t
    {
      p0,
      q0,
      H,
      D,
      Ra,
      Tdop,
      Tdopp,
      Tqop,
      Tqopp,
      Xd,
      Xdp,
      Xdpp,
      Xq,
      Xqp,
      Xqpp,
      Xl,
      S10,
      S12,
      mva,
    };

    /// Buses for a Genrou generator model
    enum class GenrouBuses : size_t
    {
      bus,
    };

    /// Signal inputs for a Genrou generator model
    enum class GenrouSignalInputs : size_t
    {
      pmech,
      efd,
    };

    /// Signal outputs for a Genrou generator model
    enum class GenrouSignalOutputs : size_t
    {
      speed,
    };

    /// Variables able to be monitored for a Genrou generator model
    enum class GenrouMonitorableVariables : size_t
    {
      ir,
      ii,
      p,
      q,
      delta,
      omega,
      speed,
    };

    /**
     * @brief Contains modeling data for a Genrou generator model.
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct GenrouData : ComponentData<real_type,
                                      index_type,
                                      GenrouParameters,
                                      GenrouBuses,
                                      GenrouSignalInputs,
                                      GenrouSignalOutputs,
                                      GenrouMonitorableVariables>
    {
      /// Units, symbols and descriptions of the GENROU data
      static constexpr typename GenrouData::Manifest MANIFEST{
          .class_name  = "Genrou",
          .description = "Round-rotor synchronous machine model",
          .parameters  = {{
              {GenrouParameters::p0, "p.u.", "P_0", "Initial active power injection"},
              {GenrouParameters::q0, "p.u.", "Q_0", "Initial reactive power injection"},
              {GenrouParameters::H, "s", "H", "Rotor inertia"},
              {GenrouParameters::D, "p.u.", "D", "Damping coefficient"},
              {GenrouParameters::Ra, "p.u.", "R_\\mathrm{a}", "Winding resistance"},
              {GenrouParameters::Tdop, "s", "T'_{d0}", "Open circuit direct axis transient time constant"},
              {GenrouParameters::Tdopp, "s", "T''_{d0}", "Open circuit direct axis sub-transient time constant"},
              {GenrouParameters::Tqop, "s", "T'_{q0}", "Open circuit quadrature axis transient time constant"},
              {GenrouParameters::Tqopp, "s", "T''_{q0}", "Open circuit quadrature axis sub-transient time constant"},
              {GenrouParameters::Xd, "p.u.", "X_d", "Direct axis synchronous reactance"},
              {GenrouParameters::Xdp, "p.u.", "X'_{d}", "Direct axis transient reactance"},
              {GenrouParameters::Xdpp, "p.u.", "X''_{d}", "Direct axis sub-transient reactance"},
              {GenrouParameters::Xq, "p.u.", "X_q", "Quadrature axis synchronous reactance"},
              {GenrouParameters::Xqp, "p.u.", "X'_{q}", "Quadrature axis transient reactance"},
              {GenrouParameters::Xqpp, "p.u.", "X''_{q}", "Quadrature axis sub-transient reactance"},
              {GenrouParameters::Xl, "p.u.", "X_{\\ell}", "Stator leakage reactance"},
              {GenrouParameters::S10, "p.u.", "S_{10}", "Saturation factor at 1.0 pu flux"},
              {GenrouParameters::S12, "p.u.", "S_{12}", "Saturation factor at 1.2 pu flux"},
              {GenrouParameters::mva, "MVA", "S^\\mathrm{base}", "Machine power base"},
          }},
          .buses       = {{
              {GenrouBuses::bus, "Terminal bus"},
          }},
          .inputs      = {{
              {GenrouSignalInputs::pmech, "Mechanical-power input"},
              {GenrouSignalInputs::efd, "Field-voltage input"},
          }},
          .outputs     = {{
              {GenrouSignalOutputs::speed, "Machine speed-deviation output"},
          }},
          .monitors    = {{
              {GenrouMonitorableVariables::ir, "p.u.", "Terminal current, real component, system base"},
              {GenrouMonitorableVariables::ii, "p.u.", "Terminal current, imaginary component, system base"},
              {GenrouMonitorableVariables::p, "p.u.", "Active power, system base"},
              {GenrouMonitorableVariables::q, "p.u.", "Reactive power, system base"},
              {GenrouMonitorableVariables::delta, "rad", "Machine internal rotor angle"},
              {GenrouMonitorableVariables::omega, "p.u.", "Machine speed deviation"},
              {GenrouMonitorableVariables::speed, "p.u.", "Per-unit machine speed"},
          }},
      };

      static_assert(MANIFEST.complete());
    };
  } // namespace PhasorDynamics
} // namespace GridKit

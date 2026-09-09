/**
 * @file FilterData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the three-phase LCL filter
 *
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Initial parameters for a filter
    enum class FilterParameters
    {
      Rs, ///< \f$\mathbf{R}_\mathrm{s}\f$ Converter-side resistance [ohm]
      Ls, ///< \f$\mathbf{L}_\mathrm{s}\f$ Converter-side inductance [H]
      C,  ///< \f$\mathbf{C}\f$ Shunt capacitance [F]
      Rg, ///< \f$\mathbf{R}_g\f$ Grid-side resistance [ohm]
      Lg, ///< \f$\mathbf{L}_g\f$ Grid-side inductance [H]
    };

    /// Inputs supported by a filter
    enum class FilterInputs : size_t
    {
      va, ///< Terminal bus phase-a voltage [V]
      vb, ///< Terminal bus phase-b voltage [V]
      vc, ///< Terminal bus phase-c voltage [V]
      ea, ///< Converter phase-a voltage [V]
      eb, ///< Converter phase-b voltage [V]
      ec, ///< Converter phase-c voltage [V]
      SIZE,
    };

    /// Outputs supported by a filter
    enum class FilterOutputs : size_t
    {
      ia,  ///< Converter-side phase-a current [A]
      ib,  ///< Converter-side phase-b current [A]
      ic,  ///< Converter-side phase-c current [A]
      voa, ///< Capacitor phase-a voltage [V]
      vob, ///< Capacitor phase-b voltage [V]
      voc, ///< Capacitor phase-c voltage [V]
      iga, ///< Phase-a current injected into the terminal bus [A]
      igb, ///< Phase-b current injected into the terminal bus [A]
      igc, ///< Phase-c current injected into the terminal bus [A]
      SIZE,
    };

    /// Variables able to be monitored for a filter
    enum class FilterMonitorableVariables
    {
      ia,
      ib,
      ic,
      voa,
      vob,
      voc,
      iga,
      igb,
      igc,
    };

    template <typename real_type, typename index_type>
    struct FilterData : public ComponentData<real_type,
                                             index_type,
                                             FilterParameters,
                                             FilterInputs,
                                             FilterOutputs,
                                             FilterMonitorableVariables>
    {
      FilterData() = default;

      using Parameters           = FilterParameters;
      using Inputs               = FilterInputs;
      using Outputs              = FilterOutputs;
      using MonitorableVariables = FilterMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit

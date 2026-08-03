/**
 * @file FrequencySweep.hpp
 *
 * @brief Frequency sweep engine shared by the EMT applications.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <cmath>

#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Parameters/Overhead.hpp>
#include <GridKit/Model/LogEvaluator.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include "FrequencyResponseData.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Application
    {
      /**
       * @brief Sweep an overhead line-parameter model over frequency.
       *
       * Poses the sweep as a DAE in s = ln(omega) and integrates it with
       * IDA, emitting one monitor sample per grid point. The model's
       * monitor configuration decides where samples go, so drivers that
       * collect samples in memory attach their own monitor sink before
       * calling. Every state is algebraic and error-controlled against
       * the model's per-variable absolute tolerances scaled by
       * ida_settings.tolerance.
       *
       * @return 0 only when the sweep reached the stop frequency and
       *         emitted exactly frequency.points monitor samples; -1 on
       *         a monitor-grid mismatch. Hard integrator failures throw
       *         SundialsException.
       */
      template <typename scalar_type, typename index_type>
      int runFrequencySweep(
          Parameters::OverheadData<scalar_type, index_type>& data,
          const FrequencyGrid&                               frequency,
          const IdaSettings&                                 ida_settings)
      {
        using namespace AnalysisManager::Sundials;

        constexpr scalar_type pi = GridKit::EMT::Constants::pi<scalar_type>();
        const scalar_type     omega_start =
            2.0 * pi * static_cast<scalar_type>(frequency.start);
        const scalar_type omega_stop =
            2.0 * pi * static_cast<scalar_type>(frequency.stop);

        Parameters::Overhead<scalar_type, index_type>         model(data);
        GridKit::Model::LogEvaluator<scalar_type, index_type> log_model(
            model, omega_start);
        log_model.allocate();

        const scalar_type s_start = std::log(omega_start);
        const scalar_type s_stop  = std::log(omega_stop);

        Ida<scalar_type, index_type> ida(&log_model);
        ida.setTolerance(static_cast<scalar_type>(ida_settings.tolerance));
        ida.setMaxSteps(static_cast<index_type>(ida_settings.max_steps));
        ida.configureSimulation();
        ida.initializeSimulation(s_start, true);

        log_model.printMonitoredVariables();
        const scalar_type dt_monitor =
            (s_stop - s_start) / static_cast<scalar_type>(frequency.points - 1);

        size_t    samples = 1;
        const int retval  = ida.runSimulation(s_stop,
                                             dt_monitor,
                                             [&samples](scalar_type)
                                             { ++samples; });
        log_model.stopMonitor();

        if (retval != 0)
        {
          return retval;
        }
        return samples == frequency.points ? 0 : -1;
      }
    } // namespace Application
  } // namespace EMT
} // namespace GridKit

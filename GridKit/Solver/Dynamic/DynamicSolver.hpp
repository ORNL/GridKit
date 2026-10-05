
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include <sundials/sundials_stepper.h>

#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{
  template <class ScalarT, typename IdxT>
  class DynamicSolver
  {
  public:
    using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

    DynamicSolver(GridKit::Model::Evaluator<ScalarT, IdxT>* model)
      : model_(model)
    {
    }

    virtual ~DynamicSolver()
    {
    }

    GridKit::Model::Evaluator<ScalarT, IdxT>* getModel()
    {
      return model_;
    }

    void setTolerance(ScalarT rel_tol)
    {
      setTolerance(rel_tol, 0);
    }

    virtual void setTolerance(ScalarT rel_tol, ScalarT abs_tol_override) = 0;
    virtual void setMaxSteps(IdxT msa)                                   = 0;

    /// This integrator as a SUNStepper over its model's state, for SUNDIALS
    /// partitioned methods. The caller owns the stepper.
    virtual SUNStepper createSUNStepper() = 0;

    /**
     * @brief Make the model's algebraic variables and derivatives consistent
     * at t, keeping its differential variables fixed.
     *
     * tout sets the scale of the step that follows.
     */
    virtual int computeConsistentState(RealT t, RealT tout) = 0;

  protected:
    GridKit::Model::Evaluator<ScalarT, IdxT>* model_;
  };

  /**
   * @brief Number of output times after t0 up to tf.
   *
   * When `dt_monitor` is nonpositive, only the final time is targeted. When
   * the final interval is epsilon-sized, it is folded into the previous
   * output.
   */
  template <typename RealT>
  int monitorStepCount(RealT t0, RealT tf, RealT dt_monitor)
  {
    if (dt_monitor <= 0.0)
    {
      return 1;
    }

    const RealT n_est   = (tf - t0) / dt_monitor;
    const RealT epsilon = std::numeric_limits<RealT>::epsilon()
                          * std::max({std::abs(t0), std::abs(tf), RealT(1.0)})
                          / dt_monitor;
    return static_cast<int>(std::ceil(n_est - epsilon));
  }

  /**
   * @brief Output time for a one-based step.
   *
   * The final output is pinned exactly to `tf` to avoid roundoff in repeated
   * time-step arithmetic.
   */
  template <typename RealT>
  RealT monitorTime(RealT t0, RealT tf, RealT dt_monitor, int step, int nsteps)
  {
    if (step == nsteps)
    {
      return tf;
    }
    return std::fma(static_cast<RealT>(step), dt_monitor, t0);
  }

} // namespace AnalysisManager

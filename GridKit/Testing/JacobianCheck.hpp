/**
 * @file JacobianCheck.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Jacobian of a model against central differences of its residual.
 */

#pragma once

#include <cmath>
#include <limits>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Central difference of the residual along y_j, with y'_j moving alpha times as fast, or none if an evaluation fails
    template <typename ScalarT, typename IdxT>
    std::vector<ScalarT> residualDifference(Model::Evaluator<ScalarT, IdxT>& model, IdxT j, ScalarT alpha)
    {
      ScalarT*      y   = model.y().getData();
      ScalarT*      yp  = model.yp().getData();
      const ScalarT y0  = y[j];
      const ScalarT yp0 = yp[j];
      const ScalarT h   = std::cbrt(std::numeric_limits<ScalarT>::epsilon()) * (ONE<ScalarT> + std::abs(y0));

      std::vector<ScalarT> difference(model.size(), ZERO<ScalarT>);
      bool                 evaluated = true;
      for (const ScalarT sign : {ONE<ScalarT>, MINUS_ONE<ScalarT>})
      {
        y[j]      = y0 + sign * h;
        yp[j]     = yp0 + sign * alpha * h;
        evaluated = evaluated && model.evaluateResidual() == 0;

        const ScalarT* f = model.getResidual().getData();
        for (IdxT i = 0; i < model.size(); ++i)
        {
          difference[i] += sign * f[i] / (TWO<ScalarT> * h);
        }
      }

      y[j]  = y0;
      yp[j] = yp0;
      if (!evaluated)
      {
        difference.clear();
      }
      return difference;
    }

    /**
     * @brief The CSR Jacobian of `model` at time `t` and its current values is
     * dF/dy + alpha dF/dy', and stores every nonzero difference
     */
    template <typename ScalarT, typename IdxT>
    bool jacobianMatchesResidual(Model::Evaluator<ScalarT, IdxT>& model, ScalarT t, ScalarT alpha, ScalarT tol)
    {
      model.updateTime(t, alpha);
      if (model.evaluateJacobian() != 0 || model.getCsrJacobian() == nullptr)
      {
        return false;
      }
      const auto jacobian = MapFromCsr(model.getCsrJacobian());

      for (IdxT j = 0; j < model.size(); ++j)
      {
        const auto difference = residualDifference(model, j, alpha);
        if (difference.size() != model.size())
        {
          return false;
        }

        for (IdxT i = 0; i < model.size(); ++i)
        {
          ScalarT stored = ZERO<ScalarT>;
          if (jacobian[i].contains(j))
          {
            stored = jacobian[i].at(j);
          }
          if (!isEqual(stored, difference[i], tol))
          {
            return false;
          }
        }
      }
      return true;
    }
  } // namespace Testing
} // namespace GridKit

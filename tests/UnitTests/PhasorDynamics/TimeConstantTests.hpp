#pragma once

#include <cmath>
#include <iostream>
#include <limits>

#include <GridKit/Testing/TestHelpers.hpp>

namespace GridKit::Testing
{
  // Probe the derivative coefficient independently of the model's right-hand side.
  template <class ModelT>
  bool implicitTimeConstant(ModelT& model, size_t row, double time_constant, double sign = -1.0)
  {
    model.allocate();
    model.tagDifferentiable();
    model.y().setToConst(0.25);
    model.yp().setToConst(0.0);
    model.evaluateResidual();
    const double initial = model.getResidual().getData()[row];

    model.yp().getData()[row] = 1.0;
    model.yp().setDataUpdated();
    model.evaluateResidual();
    const double perturbed = model.getResidual().getData()[row];
    const double tolerance = 1.0e-12 + 4.0 * std::numeric_limits<double>::epsilon() * (std::abs(initial) + std::abs(perturbed));

    const bool success = model.tag()[row] == (time_constant != 0.0)
                         && std::isfinite(initial) && std::isfinite(perturbed)
                         && isEqual(perturbed - initial, sign * time_constant, tolerance);
    if (!success)
    {
      std::cout << "Implicit time constant at row " << row << ", T = " << time_constant
                << ": tag = " << model.tag()[row]
                << ", derivative coefficient = " << perturbed - initial << '\n';
    }
    return success;
  }
} // namespace GridKit::Testing

/**
 * @file VectorFitting.cpp
 *
 * @brief Relaxed vector fitting in the real-augmented formulation.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "VectorFitting.hpp"

namespace GridKit
{
  namespace Optimization
  {
    // Implementation pending; see README.md for the algorithm stages:
    //
    //   Stage A: pole relocation - equality-constrained convex QP per pass
    //            (relaxation condition as an explicit linear constraint),
    //            relocated poles from the real-form companion eigenvalues.
    //   Stage B: coefficient identification - convex QP with optional DC,
    //            symmetry, and sample-passivity constraints.
    //   Stage C: optional refinement - NLP over pole parameters and
    //            coefficients with stability inequality constraints.
    //
    // Order search and deterministic restarts wrap the stages.
    //
    // template class VectorFitting<double, long int>;
    // template class VectorFitting<double, size_t>;
  } // namespace Optimization
} // namespace GridKit

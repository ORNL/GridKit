/**
 * @file Tgov1Enzyme.cpp
 *
 */

#include "Tgov1Impl.hpp"
#include <AutomaticDifferentiation/Enzyme/SparseWrapper.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Jacobian evaluation experimental
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::evaluateJacobian()
      {
        std::cout << "Jacobian evaluation is experimental!" << std::endl;
        GridKit::Enzyme::Sparse::EnzymeSparseModelJacobian<Tgov1<ScalarT, IdxT>, ScalarT, IdxT>(this, f_.size(), y_.data(), yp_.data(), J_);

        return 0;
      }

      // Available template instantiations
      template class Tgov1<double, long int>;
      template class Tgov1<double, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

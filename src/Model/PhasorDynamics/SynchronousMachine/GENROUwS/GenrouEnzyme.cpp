
#include "GenrouImpl.hpp"
#include <AutomaticDifferentiation/Enzyme/SparseWrapper.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Genrou..." << std::endl;
      std::cout << "Jacobian evaluation is experimental!" << std::endl;

      GridKit::Enzyme::Sparse::ModelJacobian<Genrou<ScalarT, IdxT>, ScalarT, IdxT>(this, f_.size(), y_.data(), yp_.data(), J_);

      return 0;
    }

    // Available template instantiations
    template class Genrou<double, long int>;
    template class Genrou<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit

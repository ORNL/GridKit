
#include "LoadImpl.hpp"
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
    int Load<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Load..." << std::endl;
      std::cout << "Jacobian evaluation is experimental!" << std::endl;

      GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> jac(2, 2);
      std::vector<ScalarT> y(2);
      y[0] = Vr();
      y[1] = Vi();
      GridKit::Enzyme::EnzymeSparseModelJacobian<Load<ScalarT, IdxT>, ScalarT, IdxT>(this, y.size(), y.data(), J_);

      return 0;
    }

    // Available template instantiations
    template class Load<double, long int>;
    template class Load<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit

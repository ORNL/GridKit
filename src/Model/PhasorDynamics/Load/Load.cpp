
#include "LoadImpl.hpp"
#include <AutomaticDifferentiation/Enzyme/Wrapper.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian contribution computed locally
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateJacobianLocally(const std::vector<ScalarT> x, const std::vector<ScalarT> y, GridKit::LinearAlgebra::DenseMatrix<ScalarT, IdxT>& dy_dx)
    {
      std::vector<ScalarT> v(x.size());
      std::vector<ScalarT> d_y(y.size());
      for (IdxT idy = 0; idy < y.size(); ++idy)
      {
        // Elementary vector for Jacobian-vector product
        for (IdxT idx = 0; idx < x.size(); ++idx)
        {
          v[idx] = 0.0;
        }
        v[idy] = 1.0;

        // Autodiff
        __enzyme_fwddiff<Load<ScalarT, IdxT>, ScalarT>(
            (void*) residual_wrapper<Load<ScalarT, IdxT>, ScalarT>,
            enzyme_const,
            this,
            enzyme_dup,
            x,
            v,
            enzyme_dupnoneed,
            y,
            &d_y);

        // Store result
        for (IdxT idx = 0; idx < x.size(); ++idx)
        {
          dy_dx.setValue(idx, idy, d_y[idx]);
        }
      }

      return 0;
    }

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

      GridKit::LinearAlgebra::DenseMatrix<ScalarT, IdxT> jac(2, 2);
      std::vector<ScalarT>                               y(2);
      std::vector<ScalarT>                               f(2);
      y[0] = Vr();
      y[1] = Vi();
      evaluateJacobianLocally(y, f, jac);

      J_ = *(jac.getValuesCOO()); ///< Forced setting of Jacobian
                                  ///< Careful here, because J_ was not set in the constructor

      return 0;
    }

    // Available template instantiations
    template class Load<double, long int>;
    template class Load<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit

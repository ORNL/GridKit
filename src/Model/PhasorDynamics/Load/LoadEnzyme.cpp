
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

      std::vector<ScalarT> y(2);
      std::vector<ScalarT> yp(2);
      std::vector<ScalarT> f(2);
      y[0] = Vr();
      y[1] = Vi();
      /// Setting J_ via Enzyme works, though J_ had not been initialized within the model.
      /// This is because the COO_Matrix class is very permissive.
      /// Having the currents as model variables and allocating J_ accordingly will be helpful.
      GridKit::Enzyme::Sparse::EnzymeSparseModelJacobian<Load<ScalarT, IdxT>, ScalarT, IdxT>(this, f.size(), y.data(), yp.data(), J_);

      return 0;
    }

    // Available template instantiations
    template class Load<double, long int>;
    template class Load<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit



#include "Inductor.hpp"

#include <cmath>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a inductor
   *
   * Calls default ModelEvaluatorImpl constructor.
   */

  template <class ScalarT, typename IdxT>
  Inductor<ScalarT, IdxT>::Inductor(IdxT id, RealT L, NodeT* node1, NodeT* node2)
    : L_(L), node1_(node1), node2_(node2)
  {
    assert(node1_->size() == 1);
    assert(node2_->size() == 1);
    size_           = 3;
    n_intern_       = 1;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
    nnz_            = 5;
  }

  template <class ScalarT, typename IdxT>
  Inductor<ScalarT, IdxT>::~Inductor()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Compute the absolute tolerance for each variable in the model
   *
   * @param rel_tol The relative tolerance which can be used to pick the
   *        absolute tolerance.
   * @tparam ScalarT Scalar data type
   * @tparam IdxT Index data type
   * @return int 0 if successful, non-zero otherwise.
   *
   * This represents a "noise" level close to zero for which pure relative
   * error cannot be used.
   */
  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Compute the resisdual of the component
   *
   */
  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::evaluateInternalResidual()
  {
    auto* y = y_.getData();

    f_int_[0] = -L_ * yp_int_[0] + y[1] - y[0];
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::evaluateExternalResidual()
  {
    auto* f = f_.getData();

    // input
    f[0] = -y_int_[0];
    // output
    f[1] = y_int_[0];
    return 0;
  }

  /**
   * @brief Evaluate the jacobian of the component
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    // Create dF/dy
    std::vector<IdxT>  rcord{0, 1, 2, 2, 2};
    std::vector<IdxT>  ccord{2, 2, 0, 1, 2};
    std::vector<RealT> vals{-1.0, 1.0, -1.0, 1.0, -L_ * alpha_};
    this->setJacValues(rcord, ccord, vals);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::allocate()
  {
    CircuitComponent<ScalarT, IdxT>::allocate();

    this->setExternalConnectionNodes(0, node1_->getNodeConnection(0));
    this->setExternalConnectionNodes(1, node2_->getNodeConnection(0));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Inductor<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class Inductor<double, long int>;
  template class Inductor<double, size_t>;
  template class Inductor<DependencyTracking::Variable, long int>;
  template class Inductor<DependencyTracking::Variable, size_t>;

} // namespace GridKit

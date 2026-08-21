

#include "Resistor.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a resistor model
   *
   * Calls default ModelEvaluatorImpl constructor.
   */

  template <class ScalarT, typename IdxT>
  Resistor<ScalarT, IdxT>::Resistor(IdxT id, RealT R, NodeT* node1, NodeT* node2)
    : R_(R), node1_(node1), node2_(node2)
  {
    assert(node1_->size() == 1);
    assert(node2_->size() == 1);
    size_           = 2;
    n_intern_       = 0;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
    nnz_            = 4;
  }

  template <class ScalarT, typename IdxT>
  Resistor<ScalarT, IdxT>::~Resistor()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::tagDifferentiable()
  {
    // No internal variables, so these valeus do not matter
    tag_.resize(size_, true);
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
  int Resistor<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Computes the resistors resisdual
   *
   */
  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::evaluateInternalResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::evaluateExternalResidual()
  {
    // input
    *f_ext_[0] += (*y_ext_[0] - *y_ext_[1]) / R_;
    // ouput
    *f_ext_[1] += (*y_ext_[1] - *y_ext_[0]) / R_;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    // Create dF/dy
    // does compiler make constant???
    std::vector<IdxT>  rcord{0, 0, 1, 1};
    std::vector<IdxT>  ccord{0, 1, 0, 1};
    std::vector<RealT> vals{1.0 / R_, -1.0 / R_, -1.0 / R_, 1.0 / R_};
    this->setJacValues(rcord, ccord, vals);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::allocate()
  {
    CircuitComponent<ScalarT, IdxT>::allocate();

    this->setExternalConnectionNodes(0, node1_->getNodeConnection(0));
    this->setExternalConnectionNodes(1, node2_->getNodeConnection(0));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class Resistor<double, long int>;
  template class Resistor<double, size_t>;
  template class Resistor<DependencyTracking::Variable, long int>;
  template class Resistor<DependencyTracking::Variable, size_t>;

} // namespace GridKit

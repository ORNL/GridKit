

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
  Resistor<ScalarT, IdxT>::Resistor(IdxT id, RealT R)
    : R_(R)
  {
    size_           = 2;
    n_intern_       = 0;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
  }

  template <class ScalarT, typename IdxT>
  Resistor<ScalarT, IdxT>::~Resistor()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::allocate()
  {
    y_.resize(static_cast<size_t>(size_));
    yp_.resize(static_cast<size_t>(size_));
    f_.resize(static_cast<size_t>(size_));
    tag_.resize(static_cast<size_t>(size_));
    abs_tol_.resize(static_cast<size_t>(size_));

    return 0;
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
    std::fill(abs_tol_.begin(), abs_tol_.end(), rel_tol);
    return 0;
  }

  /**
   * @brief Computes the resistors resisdual
   *
   */
  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::evaluateResidual()
  {
    // input
    f_[0] = (y_[0] - y_[1]) / R_;
    // ouput
    f_[1] = (y_[1] - y_[0]) / R_;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Resistor<ScalarT, IdxT>::evaluateJacobian()
  {

    // Create dF/dy
    // does compiler make constant???
    std::vector<IdxT>  rcord{0, 0, 1, 1};
    std::vector<IdxT>  ccord{0, 1, 0, 1};
    std::vector<RealT> vals{1.0 / R_, -1.0 / R_, -1.0 / R_, 1.0 / R_};
    jac_.setValues(rcord, ccord, vals);

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

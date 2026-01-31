
#include "BusPV.hpp"

#include <cmath>
#include <iostream>

namespace GridKit
{

  /*!
   * @brief Constructor for a PV bus
   *
   * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
   * - Number of equations = 1 (size_)
   * - Number of variables = 1 (size_)
   * - Number of quadratures = 0
   * - Number of optimization parameters = 0
   */
  template <class ScalarT, typename IdxT>
  BusPV<ScalarT, IdxT>::BusPV()
    : BaseBus<ScalarT, IdxT>(0), V_(0.0), theta0_(0.0)
  {
    // std::cout << "Create BusPV..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 1;
  }

  /*!
   * @brief Constructor for a PV bus
   *
   * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
   * - Number of equations = 1 (size_)
   * - Number of variables = 1 (size_)
   * - Number of quadratures = 0
   * - Number of optimization parameters = 0
   */
  template <class ScalarT, typename IdxT>
  BusPV<ScalarT, IdxT>::BusPV(ScalarT V, ScalarT theta0)
    : BaseBus<ScalarT, IdxT>(0), V_(V), theta0_(theta0)
  {
    // std::cout << "Create BusPV..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 1;
  }

  template <class ScalarT, typename IdxT>
  BusPV<ScalarT, IdxT>::BusPV(BusData& data)
    : BaseBus<ScalarT, IdxT>(data.bus_i), V_(data.Vm), theta0_(data.Va)
  {
    // std::cout << "Create BusPV ..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 1;
  }

  template <class ScalarT, typename IdxT>
  BusPV<ScalarT, IdxT>::~BusPV()
  {
    // std::cout << "Destroy Gen2..." << std::endl;
  }

  /*!
   * @brief allocate method resizes local solution and residual vectors.
   */
  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::allocate()
  {
    // std::cout << "Allocate PV bus ..." << std::endl;
    f_.resize(static_cast<size_t>(size_));
    y_.resize(static_cast<size_t>(size_));
    yp_.resize(static_cast<size_t>(size_));
    tag_.resize(static_cast<size_t>(size_));
    absTol_.resize(static_cast<size_t>(size_));

    fB_.resize(static_cast<size_t>(size_));
    yB_.resize(static_cast<size_t>(size_));
    ypB_.resize(static_cast<size_t>(size_));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::tagDifferentiable()
  {
    tag_[0] = false;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::setAbsoluteTolerance()
  {
    return 0;
  }

  /*!
   * @brief initialize method sets bus variables to stored initial values.
   */
  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::initialize()
  {
    // std::cout << "Initialize BusPV..." << std::endl;
    theta() = theta0_;
    yp_[0]  = 0.0;

    return 0;
  }

  /*!
   * @brief PV bus does not compute residuals, so here we just reset residual values.
   *
   * @warning This implementation assumes bus residuals are always evaluated
   * _before_ component model residuals.
   *
   */
  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::evaluateResidual()
  {
    // std::cout << "Evaluating residual of a PV bus ...\n";
    P() = 0.0; // <-- Residual P
    Q() = 0.0; // <-- Output Qg, the reactive power generator needs to supply

    return 0;
  }

  /*!
   * @brief initialize method sets bus variables to stored initial values.
   */
  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::initializeAdjoint()
  {
    // std::cout << "Initialize BusPV..." << std::endl;
    yB_[0]  = 0.0;
    ypB_[0] = 0.0;

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    fB_[0] = 0.0;

    return 0;
  }

  // Available template instantiations
  template class BusPV<double, long int>;
  template class BusPV<double, size_t>;

} // namespace GridKit

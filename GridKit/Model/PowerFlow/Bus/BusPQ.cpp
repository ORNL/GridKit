
#include "BusPQ.hpp"

#include <cmath>
#include <iostream>

namespace GridKit
{

  /*!
   * @brief Constructor for a PQ bus
   *
   * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
   * - Number of equations = 2 (size_)
   * - Number of variables = 2 (size_)
   * - Number of quadratures = 0
   * - Number of optimization parameters = 0
   */
  template <typename scalar_type, typename index_type>
  BusPQ<scalar_type, index_type>::BusPQ()
    : BaseBus<ScalarT, IdxT>(0), V0_(0.0), theta0_(0.0)
  {
    // std::cout << "Create BusPQ..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 2;
  }

  /*!
   * @brief BusPQ constructor.
   *
   * This constructor sets initial values for voltage and phase angle.
   *
   * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
   * - Number of equations = 2 (size_)
   * - Number of variables = 2 (size_)
   * - Number of quadratures = 0
   * - Number of optimization parameters = 0
   */
  template <typename scalar_type, typename index_type>
  BusPQ<scalar_type, index_type>::BusPQ(ScalarT V, ScalarT theta)
    : BaseBus<ScalarT, IdxT>(0), V0_(V), theta0_(theta)
  {
    // std::cout << "Create BusPQ..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 2;
  }

  template <typename scalar_type, typename index_type>
  BusPQ<scalar_type, index_type>::BusPQ(BusData& data)
    : BaseBus<ScalarT, IdxT>(data.bus_i), V0_(data.Vm), theta0_(data.Va)
  {
    // std::cout << "Create BusPQ..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 2;
  }

  template <typename scalar_type, typename index_type>
  BusPQ<scalar_type, index_type>::~BusPQ()
  {
    // std::cout << "Destroy PQ bus ..." << std::endl;
  }

  /*!
   * @brief allocate method resizes local solution and residual vectors.
   */
  template <typename scalar_type, typename index_type>
  int BusPQ<scalar_type, index_type>::allocate()
  {
    // std::cout << "Allocate PQ bus ..." << std::endl;
    this->allocateVectors(size_);
    tag_.resize(static_cast<size_t>(size_));

    fB_.resize(size_);
    yB_.resize(size_);
    ypB_.resize(size_);

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int BusPQ<scalar_type, index_type>::tagDifferentiable()
  {
    tag_[0] = false;
    tag_[1] = false;
    return 0;
  }

  /**
   * @brief Compute the absolute tolerance for each variable in the model
   *
   * @param rel_tol The relative tolerance which can be used to pick the
   *        absolute tolerance.
   * @return int 0 if successful, non-zero otherwise.
   *
   * This represents a "noise" level close to zero for which pure relative
   * error cannot be used.
   */
  template <typename scalar_type, typename index_type>
  int BusPQ<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /*!
   * @brief initialize method sets bus variables to stored initial values.
   */
  template <typename scalar_type, typename index_type>
  int BusPQ<scalar_type, index_type>::initialize()
  {
    // std::cout << "Initialize BusPQ..." << std::endl;
    auto* y  = y_.getData();
    auto* yp = yp_.getData();
    y[0]     = V0_;
    y[1]     = theta0_;
    yp[0]    = 0.0;
    yp[1]    = 0.0;

    y_.setDataUpdated();
    yp_.setDataUpdated();

    return 0;
  }

  /*!
   * @brief PQ bus does not compute residuals, so here we just reset residual values.
   *
   * @warning This implementation assumes bus residuals are always evaluated
   * _before_ component model residuals.
   *
   */
  template <typename scalar_type, typename index_type>
  int BusPQ<scalar_type, index_type>::evaluateResidual()
  {
    // std::cout << "Evaluating residual of a PQ bus ...\n";
    auto* f = f_.getData();
    f[0]    = 0.0;
    f[1]    = 0.0;
    f_.setDataUpdated();
    return 0;
  }

  /*!
   * @brief initialize method sets bus variables to stored initial values.
   */
  template <typename scalar_type, typename index_type>
  int BusPQ<scalar_type, index_type>::initializeAdjoint()
  {
    // std::cout << "Initialize BusPQ..." << std::endl;
    auto* yB  = yB_.getData();
    auto* ypB = ypB_.getData();
    yB[0]     = 0.0;
    yB[1]     = 0.0;
    ypB[0]    = 0.0;
    ypB[1]    = 0.0;

    yB_.setDataUpdated();
    ypB_.setDataUpdated();

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int BusPQ<scalar_type, index_type>::evaluateAdjointResidual()
  {
    auto* fB = fB_.getData();
    fB[0]    = 0.0;
    fB[1]    = 0.0;

    fB_.setDataUpdated();

    return 0;
  }

  // Available template instantiations
  template class BusPQ<double, long int>;
  template class BusPQ<double, size_t>;

} // namespace GridKit

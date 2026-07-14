
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
    using VectorT             = typename ModelEvaluatorImpl<ScalarT, IdxT>::VectorT;
    auto allocate_host_vector = [](VectorT& vector, IdxT n)
    {
      vector.resize(n);
      vector.allocate(memory::HOST);
      vector.setToZero(memory::HOST);
    };

    // std::cout << "Allocate PV bus ..." << std::endl;
    this->allocateVectors(size_);
    tag_.resize(static_cast<size_t>(size_));

    allocate_host_vector(fB_, size_);
    allocate_host_vector(yB_, size_);
    allocate_host_vector(ypB_, size_);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::tagDifferentiable()
  {
    tag_[0] = false;
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
  int BusPV<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /*!
   * @brief initialize method sets bus variables to stored initial values.
   */
  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::initialize()
  {
    // std::cout << "Initialize BusPV..." << std::endl;
    theta()  = theta0_;
    auto* yp = yp_.getData();
    yp[0]    = 0.0;

    y_.setDataUpdated();
    yp_.setDataUpdated();

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

    f_.setDataUpdated();

    return 0;
  }

  /*!
   * @brief initialize method sets bus variables to stored initial values.
   */
  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::initializeAdjoint()
  {
    // std::cout << "Initialize BusPV..." << std::endl;
    auto* yB  = yB_.getData();
    auto* ypB = ypB_.getData();
    yB[0]     = 0.0;
    ypB[0]    = 0.0;

    yB_.setDataUpdated();
    ypB_.setDataUpdated();

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPV<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    auto* fB = fB_.getData();
    fB[0]    = 0.0;

    fB_.setDataUpdated();

    return 0;
  }

  // Available template instantiations
  template class BusPV<double, long int>;
  template class BusPV<double, size_t>;

} // namespace GridKit


#include "BusSlack.hpp"

#include <cmath>
#include <iostream>

namespace GridKit
{

  /*!
   * @brief Constructor for a slack bus
   *
   * Arguments passed to ModelEvaluatorImpl:
   * - Number of equations = 0 (size_)
   * - Number of variables = 0 (size_)
   * - Number of quadratures = 0
   * - Number of optimization parameters = 0
   */
  template <class ScalarT, typename IdxT>
  BusSlack<ScalarT, IdxT>::BusSlack()
    : BaseBus<ScalarT, IdxT>(0), V_(0.0), theta_(0.0), P_(0.0), Q_(0.0), PB_(0.0), QB_(0.0)
  {
    // std::cout << "Create BusSlack..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 0;
  }

  /*!
   * @brief BusSlack constructor.
   *
   * Arguments passed to ModelEvaluatorImpl:
   * - Number of equations = 0 (size_)
   * - Number of variables = 0 (size_)
   * - Number of quadratures = 0
   * - Number of optimization parameters = 0
   */
  template <class ScalarT, typename IdxT>
  BusSlack<ScalarT, IdxT>::BusSlack(ScalarT V, ScalarT theta)
    : BaseBus<ScalarT, IdxT>(0), V_(V), theta_(theta), P_(0.0), Q_(0.0), PB_(0.0), QB_(0.0)
  {
    // std::cout << "Create BusSlack..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;
    P()   = 0.0;
    Q()   = 0.0;
    size_ = 0;
  }

  template <class ScalarT, typename IdxT>
  BusSlack<ScalarT, IdxT>::BusSlack(BusData& data)
    : BaseBus<ScalarT, IdxT>(data.bus_i), V_(data.Vm), theta_(data.Va)
  {
    // std::cout << "Create BusSlack..." << std::endl;
    // std::cout << "Number of equations is " << size_ << std::endl;
    P()   = 0.0;
    Q()   = 0.0;
    size_ = 0;
  }

  template <class ScalarT, typename IdxT>
  BusSlack<ScalarT, IdxT>::~BusSlack()
  {
  }

  template <class ScalarT, typename IdxT>
  int BusSlack<ScalarT, IdxT>::evaluateResidual()
  {
    // std::cout << "Evaluating residual of a slack bus ...\n";
    P() = 0.0;
    Q() = 0.0;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusSlack<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    PB() = 0.0;
    QB() = 0.0;
    return 0;
  }

  // Available template instantiations
  template class BusSlack<double, long int>;
  template class BusSlack<double, size_t>;

} // namespace GridKit

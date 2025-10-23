
#include "MiniGrid.hpp"

#include <cmath>
#include <iostream>
#include <vector>

#include <GridKit/Model/PowerFlow/Bus/BaseBus.hpp>

namespace GridKit
{

  /*!
   * @brief Constructor for a constant load model
   *
   * Calls default ModelEvaluatorImpl constructor.
   */

  template <class ScalarT, typename IdxT>
  MiniGrid<ScalarT, IdxT>::MiniGrid()
    : ModelEvaluatorImpl<ScalarT, IdxT>(3, 0, 0),
      Pl2_(2.5),
      Ql2_(-0.8),
      Pg3_(2.0),
      V1_(1.0),
      th1_(0.0),
      V3_(1.1),
      B12_(10.0),
      B13_(15.0),
      B22_(-22.0),
      B23_(12.0)
  {
    // std::cout << "Create a load model with " << size_ << " variables ...\n";
    rel_tol_ = 1e-5;
    abs_tol_ = 1e-5;
  }

  template <class ScalarT, typename IdxT>
  MiniGrid<ScalarT, IdxT>::~MiniGrid()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <class ScalarT, typename IdxT>
  int MiniGrid<ScalarT, IdxT>::allocate()
  {
    return 0;
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int MiniGrid<ScalarT, IdxT>::initialize()
  {
    th2() = 0.0; // th2
    V2()  = 1.0; // V2
    th3() = 0.0; // th3
    return 0;
  }

  /**
   * @brief Contributes to the bus residual.
   *
   * Must be connected to a PQ bus.
   */
  template <class ScalarT, typename IdxT>
  int MiniGrid<ScalarT, IdxT>::evaluateResidual()
  {
    f_[0] = -Pl2_ - V2() * (V1_ * B12_ * sin(th2() - th1_) + V3_ * B23_ * sin(th2() - th3()));
    f_[1] = -Ql2_ + V2() * (V1_ * B12_ * cos(th2() - th1_) + B22_ * V2() + V3_ * B23_ * cos(th2() - th3()));
    f_[2] = Pg3_ - V3_ * (V1_ * B13_ * sin(th3() - th1_) + V2() * B23_ * sin(th3() - th2()));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MiniGrid<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  // Available template instantiations
  template class MiniGrid<double, long int>;
  template class MiniGrid<double, size_t>;

} // namespace GridKit

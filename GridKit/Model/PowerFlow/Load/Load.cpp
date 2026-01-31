
#include "Load.hpp"

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
  Load<ScalarT, IdxT>::Load(bus_type* bus, ScalarT P, ScalarT Q)
    : P_(P),
      Q_(Q),
      busID_(0),
      bus_(bus)
  {
    // std::cout << "Create a load model with " << size_ << " variables ...\n";
    size_ = 0;
  }

  template <class ScalarT, typename IdxT>
  Load<ScalarT, IdxT>::Load(bus_type* bus, LoadData& data)
    : P_(data.Pd),
      Q_(data.Qd),
      busID_(data.bus_i),
      bus_(bus)
  {
    // std::cout << "Create a load model with " << size_ << " variables ...\n";
    size_ = 0;
  }

  template <class ScalarT, typename IdxT>
  Load<ScalarT, IdxT>::~Load()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::allocate()
  {
    return 0;
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * \brief Specify absolute tolerance.
   */
  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::setAbsoluteTolerance()
  {
    return 0;
  }

  /**
   * @brief Contributes to the bus residual.
   *
   * Must be connected to a PQ bus.
   */
  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::evaluateResidual()
  {
    // std::cout << "Evaluating load residual ...\n";
    bus_->P() -= P_;
    bus_->Q() -= Q_;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Load<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class Load<double, long int>;
  template class Load<double, size_t>;

} // namespace GridKit

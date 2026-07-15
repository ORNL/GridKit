
#include "GeneratorPV.hpp"

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
  GeneratorPV<ScalarT, IdxT>::GeneratorPV(bus_type* bus, GenData& data)
    : P_(data.Pg),
      // Q_(data.Qg),
      bus_(bus)
  {
    // std::cout << "Create a load model with " << size_ << " variables ...\n";
    size_ = 0;
  }

  template <class ScalarT, typename IdxT>
  GeneratorPV<ScalarT, IdxT>::~GeneratorPV()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::allocate()
  {
    return 0;
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Contributes to the bus residual.
   *
   * Must be connected to a PQ bus.
   */
  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::evaluateResidual()
  {
    // std::cout << "Evaluating load residual ...\n";
    bus_->P() += P_;
    // bus_->Q() += Q_;
    if (bus_->size() > 0)
    {
      bus_->getResidual().setDataUpdated();
    }
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int GeneratorPV<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class GeneratorPV<double, long int>;
  template class GeneratorPV<double, size_t>;

} // namespace GridKit


#include "GeneratorSlack.hpp"

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

  template <typename scalar_type, typename index_type>
  GeneratorSlack<scalar_type, index_type>::GeneratorSlack(BusT* bus, GenData& /* data */)
    : bus_(bus)
  {
    // std::cout << "Create a load model with " << size_ << " variables ...\n";
    size_ = 0;
  }

  template <typename scalar_type, typename index_type>
  GeneratorSlack<scalar_type, index_type>::~GeneratorSlack()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::allocate()
  {
    return 0;
  }

  /**
   * Initialization of the grid model
   */
  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Contributes to the bus residual.
   *
   * Must be connected to a PQ bus.
   */
  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::evaluateResidual()
  {
    // std::cout << "Evaluating load residual ...\n";
    // bus_->P() += P_;
    // bus_->Q() += Q_;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::evaluateJacobian()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::evaluateIntegrand()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::initializeAdjoint()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int GeneratorSlack<scalar_type, index_type>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class GeneratorSlack<double, long int>;
  template class GeneratorSlack<double, size_t>;

} // namespace GridKit

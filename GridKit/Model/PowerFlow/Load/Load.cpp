
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

  template <typename scalar_type, typename index_type>
  Load<scalar_type, index_type>::Load(BusT* bus, ScalarT P, ScalarT Q)
    : P_(P),
      Q_(Q),
      busID_(0),
      bus_(bus)
  {
    // std::cout << "Create a load model with " << size_ << " variables ...\n";
    size_ = 0;
  }

  template <typename scalar_type, typename index_type>
  Load<scalar_type, index_type>::Load(BusT* bus, LoadData& data)
    : P_(data.Pd),
      Q_(data.Qd),
      busID_(data.bus_i),
      bus_(bus)
  {
    // std::cout << "Create a load model with " << size_ << " variables ...\n";
    size_ = 0;
  }

  template <typename scalar_type, typename index_type>
  Load<scalar_type, index_type>::~Load()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::allocate()
  {
    return 0;
  }

  /**
   * Initialization of the grid model
   */
  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::tagDifferentiable()
  {
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
  int Load<scalar_type, index_type>::setAbsoluteTolerance(RealT)
  {
    return 0;
  }

  /**
   * @brief Contributes to the bus residual.
   *
   * Must be connected to a PQ bus.
   */
  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::evaluateResidual()
  {
    // std::cout << "Evaluating load residual ...\n";
    bus_->P() -= P_;
    bus_->Q() -= Q_;
    if (bus_->size() > 0)
    {
      bus_->getResidual().setDataUpdated();
    }
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::evaluateJacobian()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::evaluateIntegrand()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::initializeAdjoint()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Load<scalar_type, index_type>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class Load<double, long int>;
  template class Load<double, size_t>;

} // namespace GridKit

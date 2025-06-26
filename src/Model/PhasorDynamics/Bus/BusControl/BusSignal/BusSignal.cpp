
/**
 * @file BusSignal.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of BusSignal class.
 *
 */

#include "BusSignal.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/BusData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {

        /**
     * @brief Construct a new Bus
     *
     * @tparam ScalarT - type of scalar variables
     * @tparam IdxT    - type for vector/matrix indices
     * @param[in] data - structure with bus data
     */
    template <class ScalarT, typename IdxT>
    BusSignal<ScalarT, IdxT>::BusSignal()
      : BusSignalBase<ScalarT, IdxT>()
    {
      // std::cout << "Create Bus..." << std::endl;
      // std::cout << "Number of equations is " << size_ << std::endl;

      size_ = 0;
    }

    /**
     * @brief Construct a new Bus
     *
     * @tparam ScalarT - type of scalar variables
     * @tparam IdxT    - type for vector/matrix indices
     * @param[in] data - structure with bus data
     */
    template <class ScalarT, typename IdxT>
    BusSignal<ScalarT, IdxT>::BusSignal(const DataT& data)
      : BusSignalBase<ScalarT, IdxT>(data.bus_id)
    {
      // std::cout << "Create Bus..." << std::endl;
      // std::cout << "Number of equations is " << size_ << std::endl;

      size_ = 0;
    }

    template <class ScalarT, typename IdxT>
    BusSignal<ScalarT, IdxT>::~BusSignal()
    {
      // std::cout << "Destroy PQ bus ..." << std::endl;
    }

    /*!
     * @brief allocate method resizes local solution and residual vectors.
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::allocate()
    {
      // Temporary while we use std::vector in the code
      size_t size = static_cast<size_t>(size_);

      // Resize component model data
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      fB_.resize(size);
      yB_.resize(size);
      ypB_.resize(size);
      return 0;
    }

    /*!
     * @brief Bus variables are algebraic.
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    /*!
     * @brief initialize method sets bus variables to stored initial values.
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::initialize()
    {
      // std::cout << "Initialize Bus..." << std::endl;
      return 0;
    }

    /*!
     * @brief PQ bus does not compute residuals, so here we just reset residual values.
     *
     * @warning This implementation assumes bus residuals are always evaluated
     * _before_ component model residuals.
     *
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateResidual()
    {
      // std::cout << "Evaluating residual of a PQ bus ...\n";
      return 0;
    }

    /**
     * @brief Jacobian evaluation not implemented
     *
     * @tparam ScalarT - data type for Jacobian elements
     * @tparam IdxT    - data type for matrix indices
     * @return int - error code
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateJacobian()
    {
      return 0;
    }

    /*!
     * @brief initialize method sets bus variables to stored initial values.
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::initializeAdjoint()
    {
      // std::cout << "Initialize Bus..." << std::endl;
      return 0;
    }

    /**
     * @brief Bus only initializes adjoint residual elements to zero.
     *
     * @tparam ScalarT - data type for the integrand
     * @tparam IdxT    - data type for matrix/vector indices
     * @return int - error code
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateAdjointResidual()
    {
      return 0;
    }

    /**
     * @brief Quadrature evaluation not implemented
     *
     * @tparam ScalarT - data type for the integrand
     * @tparam IdxT    - data type for matrix/vector indices
     * @return int - error code
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateIntegrand()
    {
      return 0;
    }

    /**
     * @brief Adjoint quadrature evaluation not implemented
     *
     * @tparam ScalarT - data type for the integrand
     * @tparam IdxT    - data type for matrix/vector indices
     * @return int - error code
     */
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
      return 0;
    }

  } // namespace PhasorDynamics
} // namespace GridKit

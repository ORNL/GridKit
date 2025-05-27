/**
 * @file BusFault.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a bus fault to ground model.
 *
 * The model uses Cartesian coordinates.
 *
 */

#include "BusFault.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PowerFlow/PowerSystemData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Constructor for a pi-model branch
     *
     * Arguments passed to ModelEvaluatorImpl:
     * - Number of equations = 0
     * - Number of independent variables = 0
     * - Number of quadratures = 0
     * - Number of optimization parameters = 0
     */
    template <class ScalarT, typename IdxT>
    BusFault<ScalarT, IdxT>::BusFault(bus_type* bus)
      : bus_(bus), R_(0), X_(0.01), status_(0), busID_(0)
    {
      (void) busID_;
      size_ = 0;
    }

    /**
     * @brief Construct a new BusFault
     *
     * @tparam ScalarT - scalar type
     * @tparam IdxT    - matrix/vector index type
     * @param bus1 - pointer to bus-1
     * @param bus2 - pointer to bus-2
     * @param R - line series resistance
     * @param X - line series reactance
     * @param G - line shunt conductance
     * @param B - line shunt charging
     */
    template <class ScalarT, typename IdxT>
    BusFault<ScalarT, IdxT>::BusFault(bus_type* bus, real_type R, real_type X, int status)
      : bus_(bus), R_(R), X_(X), status_(status), busID_(0)
    {
      size_ = 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::allocate()
    {
      // std::cout << "Allocate BusFault..." << std::endl;
      return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * \brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     *
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::evaluateResidual()
    {
      if (status_)
      {
        double B = -X_ / (X_ * X_ + R_ * R_);
        double G = R_ / (X_ * X_ + R_ * R_);

        Ir() += -Vr() * G + Vi() * B;
        Ii() += -Vr() * B - Vi() * G;
      }
      return 0;
    }

    /**
     * @brief Jacobian evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for BusFault..." << std::endl;
      std::cout << "Jacobian evaluation not implemented!" << std::endl;
      return 0;
    }

    /**
     * @brief Integrand (objective) evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::evaluateIntegrand()
    {
      // std::cout << "Evaluate Integrand for BusFault..." << std::endl;
      return 0;
    }

    /**
     * @brief Adjoint initialization not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::initializeAdjoint()
    {
      // std::cout << "Initialize adjoint for BusFault..." << std::endl;
      return 0;
    }

    /**
     * @brief Adjoint residual evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::evaluateAdjointResidual()
    {
      // std::cout << "Evaluate adjoint residual for BusFault..." << std::endl;
      return 0;
    }

    /**
     * @brief Adjoint integrand (objective) evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
      // std::cout << "Evaluate adjoint Integrand for BusFault..." << std::endl;
      return 0;
    }

    // Available template instantiations
    template class BusFault<double, long int>;
    template class BusFault<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit

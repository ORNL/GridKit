/**
 * @file SynchronousMachine.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */

#include "SynchronousMachine.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <PowerSystemData.hpp>

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
    SynchronousMachine<ScalarT, IdxT>::SynchronousMachine(bus_type* bus)
        : bus_(bus), R_(0.0), X_(0.01), G_(0.0), B_(0.0)
    {
        size_ = 0;
    }

    /**
     * @brief Destroy the SynchronousMachine
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    SynchronousMachine<ScalarT, IdxT>::~SynchronousMachine()
    {
        // std::cout << "Destroy SynchronousMachine..." << std::endl;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int SynchronousMachine<ScalarT, IdxT>::allocate()
    {
        // std::cout << "Allocate SynchronousMachine..." << std::endl;
        return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <class ScalarT, typename IdxT>
    int SynchronousMachine<ScalarT, IdxT>::initialize()
    {
        return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int SynchronousMachine<ScalarT, IdxT>::tagDifferentiable()
    {
        return 0;
    }

    /**
     * \brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     *
     */
    template <class ScalarT, typename IdxT>
    int SynchronousMachine<ScalarT, IdxT>::evaluateResidual()
    {
        // std::cout << "Evaluating branch residual ...\n";
        // real_type b = -X_/(R_*R_ + X_*X_);
        // real_type g =  R_/(R_*R_ + X_*X_);

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
    int SynchronousMachine<ScalarT, IdxT>::evaluateJacobian()
    {
        std::cout << "Evaluate Jacobian for SynchronousMachine..." << std::endl;
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
    int SynchronousMachine<ScalarT, IdxT>::evaluateIntegrand()
    {
        // std::cout << "Evaluate Integrand for SynchronousMachine..." << std::endl;
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
    int SynchronousMachine<ScalarT, IdxT>::initializeAdjoint()
    {
        // std::cout << "Initialize adjoint for SynchronousMachine..." << std::endl;
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
    int SynchronousMachine<ScalarT, IdxT>::evaluateAdjointResidual()
    {
        // std::cout << "Evaluate adjoint residual for SynchronousMachine..." << std::endl;
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
    int SynchronousMachine<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
        // std::cout << "Evaluate adjoint Integrand for SynchronousMachine..." << std::endl;
        return 0;
    }

    // Available template instantiations
    template class SynchronousMachine<double, long int>;
    template class SynchronousMachine<double, size_t>;

} // namespace PhasorDynamics
} // namespace GridKit


#include "Load.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Load/LoadData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Constructor for a pi-model load
     *
     * Arguments passed to ModelEvaluatorImpl:
     * - Number of equations = 0
     * - Number of independent variables = 0
     * - Number of quadratures = 0
     * - Number of optimization parameters = 0
     */

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type* bus)
      : bus_(bus)
    {
      size_ = 0;
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type* bus,
                              real_type R,
                              real_type X)
      : bus_(bus),
        R_(R),
        X_(X)
    {
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type* bus, IdxT component_id)
      : bus_(bus)
    {
      size_         = 0;
      component_id_ = component_id;
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::~Load()
    {
      // std::cout << "Destroy Load..." << std::endl;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::allocate()
    {
      // std::cout << "Allocate Load..." << std::endl;
      return 0;
    }

    /**
     * Initialization of the load model
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Residual contribution computed locally
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateResidualLocally(const std::vector<ScalarT> y, std::vector<ScalarT>& f)
    {
      real_type b = -X_ / (R_ * R_ + X_ * X_);
      real_type g = R_ / (R_ * R_ + X_ * X_);

      f[0] = -g * y[0] + b * y[1];
      f[1] = -b * y[0] - g * y[1];

      return 0;
    }

    /**
     * \brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateResidual()
    {
      std::vector<ScalarT> y(2);
      std::vector<ScalarT> f(2);
      y[0] = Vr();
      y[1] = Vi();
      evaluateResidualLocally(y, f);
      Ir() += f[0];
      Ii() += f[1];

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
    int Load<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Load..." << std::endl;
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
    int Load<ScalarT, IdxT>::evaluateIntegrand()
    {
      // std::cout << "Evaluate Integrand for Load..." << std::endl;
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
    int Load<ScalarT, IdxT>::initializeAdjoint()
    {
      // std::cout << "Initialize adjoint for Load..." << std::endl;
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
    int Load<ScalarT, IdxT>::evaluateAdjointResidual()
    {
      // std::cout << "Evaluate adjoint residual for Load..." << std::endl;
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
    int Load<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
      // std::cout << "Evaluate adjoint Integrand for Load..." << std::endl;
      return 0;
    }

    // Available template instantiations
    template class Load<double, long int>;
    template class Load<double, size_t>;
    template class Load<DependencyTracking::Variable, long int>;
    template class Load<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit


#include <iostream>
#include <cmath>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <PowerSystemData.hpp>

#include "Branch.hpp"

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
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2)
    : bus1_(bus1),
      bus2_(bus2),
      R_(0.0),
      X_(0.01),
      G_(0.0),
      B_(0.0),
      bus1ID_(0),
      bus2ID_(0)
    {
        size_ = 0;
    }

    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1,
                                  bus_type* bus2,
                                  real_type R,
                                  real_type X,
                                  real_type G,
                                  real_type B)
    : bus1_(bus1),
      bus2_(bus2),
      R_(R),
      X_(X),
      G_(G),
      B_(B),
      bus1ID_(0),
      bus2ID_(0)
    {
    }

    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2, BranchData& data)
    : bus1_(bus1),
      bus2_(bus2),
      R_(data.r),
      X_(data.x),
      G_(0.0),
      B_(data.b),
      bus1ID_(data.fbus),
      bus2ID_(data.tbus)
    {
        size_ = 0;
    }


    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::~Branch()
    {
        //std::cout << "Destroy Branch..." << std::endl;
    }

    /*!
    * @brief allocate method computes sparsity pattern of the Jacobian.
    */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::allocate()
    {
        //std::cout << "Allocate Branch..." << std::endl;
        return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::initialize()
    {
        return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::tagDifferentiable()
    {
        return 0;
    }

    /**
     * \brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     * 
     * @todo Add and verify conductance to ground (B and G)
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateResidual()
    {
        // std::cout << "Evaluating branch residual ...\n";
        real_type b = -X_/(R_*R_ + X_*X_);
        real_type g =  R_/(R_*R_ + X_*X_);

        Ir1() += -(g + 0.5*G_)*Vr1() + (b + 0.5*B_)*Vi1() + g*Vr2() - b*Vi2();
        Ii1() += -(b + 0.5*B_)*Vr1() - (g + 0.5*G_)*Vi1() + b*Vr2() + g*Vi2();
        Ir2() +=  g*Vr1() - b*Vi1() - (g + 0.5*G_)*Vr2() + (b + 0.5*B_)*Vi2();
        Ii2() +=  b*Vr1() + g*Vi1() - (b + 0.5*B_)*Vr2() - (g + 0.5*G_)*Vi2();

        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateJacobian()
    {
        std::cout << "Evaluate Jacobian for Branch..." << std::endl;
        std::cout << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateIntegrand()
    {
        // std::cout << "Evaluate Integrand for Branch..." << std::endl;
        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::initializeAdjoint()
    {
        //std::cout << "Initialize adjoint for Branch..." << std::endl;
        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateAdjointResidual()
    {
        // std::cout << "Evaluate adjoint residual for Branch..." << std::endl;
        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
        // std::cout << "Evaluate adjoint Integrand for Branch..." << std::endl;
        return 0;
    }

    // Available template instantiations
    template class Branch<double, long int>;
    template class Branch<double, size_t>;

} //namespace PhasorDynamics
} //namespace GridKit

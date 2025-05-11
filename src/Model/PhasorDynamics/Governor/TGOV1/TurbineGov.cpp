/**
 * @file TurbineGov.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Definition of a Turbine Governor Model (IEEET1).
 *
 *
 */

#include "TurbineGov.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Constructor for Governor with default parameters
     *
     * @param unit_id Associated Generator Unit ID
     */
    template <class ScalarT, typename IdxT>
    TurbineGov<ScalarT, IdxT>::TurbineGov(machine_type* machine)
        : unit_id_(unit_id),
          R_(3.),
          Pvmin_(0.),
          Pvmax_(1.),
          T1_(7.),
          T2_(.04),
          T3_(.05),
          Dt_(.75),
    {
      // 2 Internal Variables and 7 Parameters
      size_       = 2;
      size_param_ = 7;
    }

    /*!
     * @brief Constructor for Governor
     *
     * @param unit_id Generator Unit ID
     * @param R       Droop Constant (p.u.) 
     * @param Pvmin   Minimum Valve Position
     * @param Pvmax   Maximum Valve Position
     * @param Dt      Turbine Damping Coefficient (p.u.)
     */
    template <class ScalarT, typename IdxT>
    TurbineGov<ScalarT, IdxT>::TurbineGov(
        machine_type* machine,
        real_type R,
        real_type Pvmin,
        real_type T1,
        real_type T2,
        real_type T3,
        real_type Dt
    )
      : machine_(machine),
        R_(3.),
        Pvmin_(0.),
        Pvmax_(1.),
        T1_(7.),
        T2_(.04),
        T3_(.05),
        Dt_(.75),
    {

      // 2 Internal Variables and 7 Parameters
      size_       = 2;
      size_param_ = 7;

    }

    /*!
     * @brief Destructor for Governor Model
     *
     */
    template <class ScalarT, typename IdxT>
    TurbineGov<ScalarT, IdxT>::~TurbineGov()
    {
      std::cout << "Destroy Governor ..." << std::endl;
    }


    /*!
     * @brief Allocate memory for model
     *
     *
     */
    template <class ScalarT, typename IdxT>
    int TurbineGov<ScalarT, IdxT>::allocate()
    {
      f_.resize(size_);
      y_.resize(size_);
      yp_.resize(size_);
      tag_.resize(size_);
      fB_.resize(size_);
      yB_.resize(size_);
      ypB_.resize(size_);
      return 0;
    }

    /**
     * @brief Initialization of the Governor
     *
     */
    template <class ScalarT, typename IdxT>
    int TurbineGov<ScalarT, IdxT>::initialize()
    {

      // Initialize Differential Variables:
      y_[0] = pref_;  // Valve Position
      y_[1] = pmech_; // Turbine Power

      // Initialize Differential Variables
      yp_[0] = 0.0;
      yp_[1] = 0.0;

      return 0;
    }

    /**
     * @brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int TurbineGov<ScalarT, IdxT>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[i] = true;
      }
      return 0;
    }

    /**
     * @brief Residuals of system equations
     *
     */
    template <class ScalarT, typename IdxT>
    int TurbineGov<ScalarT, IdxT>::evaluateResidual()
    {

      // External Variables
      // TODO: should be deviation, not true speed.
      ScalarT omega  = machine_->omega();
      
      // Internal Variables 
      ScalarT pv       = y_[0];
      ScalarT ptx      = y_[1];

      // Internal Derivatives
      ScalarT pv_dot   = yp_[0];
      ScalarT ptx_dot  = yp_[1];

      // Internal Differential Equations
      f_[0] = ptx_dot - pv + (ptx + T2_ * pv) / T3_;
      f_[1] = pv_dot  + (pv - (pref_ - omega) / R ) / T3_;
      
      // Output Variable
      pmech_ = (ptx + T2_ * pv) / T3_ - (Dt_ * omega);

      return 0;
    }

    /**
     * @brief Jacobian evaluation not implemented yet
     *
     * @tparam ScalarT - Scalar data type
     * @tparam IdxT    - Index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int TurbineGov<ScalarT, IdxT>::evaluateJacobian()
    {
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
    int TurbineGov<ScalarT, IdxT>::evaluateIntegrand()
    {
      std::cout << "Evaluate Integrand for TurbineGov..." << std::endl;
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
    int TurbineGov<ScalarT, IdxT>::initializeAdjoint()
    {
      std::cout << "Initialize adjoint for TurbineGov..." << std::endl;
      return 0;
    }

    /**
     * @brief Adjoint residual evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int     - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int TurbineGov<ScalarT, IdxT>::evaluateAdjointResidual()
    {
      std::cout << "Evaluate adjoint residual for TurbineGov..." << std::endl;
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
    int TurbineGov<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
      std::cout << "Evaluate adjoint Integrand for TurbineGov..." << std::endl;
      return 0;
    }

    /**
     * @brief The mechanical power output.
     * 
     * @return ScalarT - Mechanical output power value.
     */
    template <class ScalarT, typename IdxT>
    ScalarT TurbineGov<ScalarT, IdxT>::Pmech()
    {
      return pmech_;
    }

    // Available template instantiations
    template class TurbineGov<double, long int>;
    template class TurbineGov<double, size_t>;

    

  } // namespace PhasorDynamics
} // namespace GridKit

/**
 * @file TurbineGov.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Definition of a Turbine Governor Model (IEEET1).
 *
 */

#include "TurbineGov.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Governor/TGOV1/TurbineGovData.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {

    /*!
     * @brief Constructor for Governor
     *
     * @param machine Generator Object
     * @param data TGOV1 Data Object
     */
    template <class ScalarT, typename IdxT>
    TurbineGov<ScalarT, IdxT>::TurbineGov(machine_type* machine, const model_data_type& data)
      : machine_(machine),
        R_(data.R),
        Pvmin_(data.Pvmin),
        Pvmax_(data.Pvmax),
        T1_(data.T1),
        T2_(data.T2),
        T3_(data.T3),
        Dt_(data.Dt)
    {

      // 3 Internal Variables
      size_ = 3;
    }

    /*!
     * @brief Allocate memory for model
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

      // Input Variables (Parameter for now)
      pref_ = 1 / R_;

      // Differential Variables
      y_[0] = (T3_ - T2_) * R_ * pref_; // Ptx (Turbine Power )
      y_[1] = R_ * pref_;               // Pv  (Valve Position)

      // Algebraic Variables
      y_[2] = R_ * pref_; // Pmech

      // D.V. Derivative
      yp_[0] = 0.0; // Ptx
      yp_[1] = 0.0; // Pv

      return 0;
    }

    /**
     * @brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int TurbineGov<ScalarT, IdxT>::tagDifferentiable()
    {

      tag_[0] = true;  // Pv
      tag_[1] = true;  // Ptx
      tag_[2] = false; // Pmech

      return 0;
    }

    /**
     * @brief Residuals of system equations
     *
     */
    template <class ScalarT, typename IdxT>
    int TurbineGov<ScalarT, IdxT>::evaluateResidual()
    {

      // Input Variables
      ScalarT omega = machine_->speed();

      // Internal Variables
      ScalarT pv    = y_[0];
      ScalarT ptx   = y_[1];
      ScalarT pmech = y_[2];

      // Internal Derivatives
      ScalarT pv_dot  = yp_[0];
      ScalarT ptx_dot = yp_[1];

      // Internal Differential Equations
      f_[0] = ptx_dot - pv + (ptx + T2_ * pv) / T3_;
      f_[1] = pv_dot + (pv - (pref_ - omega) / R_) / T3_;

      // Internal Algebraic Equations
      f_[2] = pmech - (ptx + T2_ * pv) / T3_ - (Dt_ * omega);

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
      return y_[2];
    }

    // Available template instantiations
    template class TurbineGov<double, long int>;
    template class TurbineGov<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit

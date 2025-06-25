/**
 * @file Tgov1.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Definition of a Turbine Governor Model (IEEET1).
 *
 */

#include "Tgov1.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/BusSignal/BusSignal.hpp>
#include <Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {

      /*!
       * @brief Constructor for Governor
       *
       * @param data    TGOV1 Data Object
       */
      template <class ScalarT, typename IdxT>
      Tgov1<ScalarT, IdxT>::Tgov1(const model_data_type& data)
        : speed_signal_(nullptr),
          pmech_signal_(nullptr),
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

      template <class ScalarT, typename IdxT>
      Tgov1<ScalarT, IdxT>::Tgov1()
        : speed_signal_(nullptr),
          pmech_signal_(nullptr),
          R_(0.05),
          Pvmin_(0),
          Pvmax_(1),
          T1_(0.5),
          T2_(2.5),
          T3_(7.5),
          Dt_(0)
      {

        // 3 Internal Variables
        size_ = 3;
      }

      /*!
       * @brief Allocate memory for model
       *
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::allocate()
      {
        auto size = static_cast<size_t>(size_); // avoid compiler warnings
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        fB_.resize(size);
        yB_.resize(size);
        ypB_.resize(size);
        return 0;
      }

      /**
       * @brief Initialization of the Governor
       *
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::initialize()
      {

        // Requires machine to be initialized first
        // (Is there a better way to have dependent initialization schemes?)
        ScalarT p0 = 0;
        if (pmech_signal_)
        {
          p0 = pmech_signal_->Vr();
        }

        // Input Variables (Parameter for now)
        pref_ = R_ * p0;

        // Internal States
        y_[0] = (T3_ - T2_) * p0; // y0 - Ptx (Turbine Power )
        y_[1] = p0;               // y1 - Pv  (Valve Position)
        y_[2] = p0;               // y2 - Pm  (Mech Power)

        // D.V. Derivative
        yp_[0] = 0.0; // Ptx
        yp_[1] = 0.0; // Pv
        yp_[2] = 0.0; // Pm

        return 0;
      }

      /**
       * @brief Identify differential variables.
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::tagDifferentiable()
      {

        tag_[0] = true;  // Pv
        tag_[1] = true;  // Ptx
        tag_[2] = false; // Pmech

        return 0;
      }

      /**
       * @brief Scaled sigmoid activation function
       *
       * Temporary local implementation of smooth approximation
       * of a piecewise differential equation. Ideally this is
       * a more abstracted capability with GK.
       *
       * Algebraic approximation of transcendental sigmoid.
       */
      template <class ScalarT, typename IdxT>
      ScalarT Tgov1<ScalarT, IdxT>::sigmoid(ScalarT x)
      {
        return ((0.5 * mu_ * x) / (1.0 + std::abs(mu_ * x))) + 0.5;
      }

      /**
       * @brief Indicator function for lower valve limit violation.
       */
      template <class ScalarT, typename IdxT>
      ScalarT Tgov1<ScalarT, IdxT>::indicator_low(ScalarT x, ScalarT f)
      {
        return (this->sigmoid(Pvmin_ - x)) * (this->sigmoid(-f));
      }

      /**
       * @brief Indicator function for high valve limit violation.
       */
      template <class ScalarT, typename IdxT>
      ScalarT Tgov1<ScalarT, IdxT>::indicator_high(ScalarT x, ScalarT f)
      {
        return (this->sigmoid(x - Pvmax_)) * (this->sigmoid(f));
      }

      /**
       * @brief Net Indicator function for valve limits.
       */
      template <class ScalarT, typename IdxT>
      ScalarT Tgov1<ScalarT, IdxT>::indicator(ScalarT x, ScalarT f)
      {
        return (1 - this->indicator_low(x, f)) * (1 - this->indicator_high(x, f));
      }

      /**
       * @brief Residuals of system equations
       *
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::evaluateResidual()
      {

        // Input Variables
        ScalarT omega = 0;
        if (speed_signal_)
        {
          omega = speed_signal_->Vr();
        }

        // Read Internal Variables
        ScalarT ptx   = y_[0]; // y0 - Ptx
        ScalarT pv    = y_[1]; // y1 - Pv
        ScalarT pmech = y_[2]; // y2 - Pmech

        // Read Internal Derivatives
        ScalarT ptx_dot = yp_[0];
        ScalarT pv_dot  = yp_[1];

        // The 'pre-limit' derivative of Pv
        ScalarT f        = (-pv + (pref_ - omega) / R_) / T1_;
        ScalarT valv_ind = this->indicator(pv, f);

        // Internal Differential Equations
        f_[0] = -ptx_dot + pv - (ptx + T2_ * pv) / T3_;
        f_[1] = -pv_dot + valv_ind * f;

        // Internal Algebraic Equations
        f_[2] = -pmech + (ptx + T2_ * pv) / T3_ - (Dt_ * omega);

        // Update signal if available
        if (pmech_signal_)
        {
          pmech_signal_->Vr() = pmech;
        }

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
      int Tgov1<ScalarT, IdxT>::evaluateJacobian()
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
      int Tgov1<ScalarT, IdxT>::evaluateIntegrand()
      {
        std::cout << "Evaluate Integrand for Tgov1..." << std::endl;
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
      int Tgov1<ScalarT, IdxT>::initializeAdjoint()
      {
        std::cout << "Initialize adjoint for Tgov1..." << std::endl;
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
      int Tgov1<ScalarT, IdxT>::evaluateAdjointResidual()
      {
        std::cout << "Evaluate adjoint residual for Tgov1..." << std::endl;
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
      int Tgov1<ScalarT, IdxT>::evaluateAdjointIntegrand()
      {
        std::cout << "Evaluate adjoint Integrand for Tgov1..." << std::endl;
        return 0;
      }

      /**
       * @brief The machine speed signal setter.
       *
       * @param signal A BusSignal pointer
       */
      template <class ScalarT, typename IdxT>
      void Tgov1<ScalarT, IdxT>::set_speed_signal(bus_type* signal)
      {
        speed_signal_ = signal;
      }

      /**
       * @brief The machine pmech signal setter.
       *
       * @param signal A BusSignal pointer
       */
      template <class ScalarT, typename IdxT>
      void Tgov1<ScalarT, IdxT>::set_pmech_signal(bus_type* signal)
      {
        pmech_signal_ = signal;
      }

      // Available template instantiations
      template class Tgov1<double, long int>;
      template class Tgov1<double, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

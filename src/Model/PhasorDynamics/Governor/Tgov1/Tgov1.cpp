/**
 * @file Tgov1.cpp
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Definition of a Turbine Governor Model (IEEET1).
 *
 */

#include "Tgov1.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#include <Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       *
       */
      template <class ScalarT, typename IdxT>
      Tgov1<ScalarT, IdxT>::Tgov1(machine_type* machine, const model_data_type& data)
        : machine_(machine)
      {
        initializeParameters(data);
        size_ = 3;
      }

      template <class ScalarT, typename IdxT>
      Tgov1<ScalarT, IdxT>::Tgov1(signal_type* pmech, signal_type* omega, const model_data_type& data)
        : pmech_(pmech),
          omega_(omega)
      {
        initializeParameters(data);
        size_ = 3;
      }

      template <class ScalarT, typename IdxT>
      void Tgov1<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
      {
        if (data.parameters.contains(model_data_type::Parameters::R))
        {
          R_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::R));
        }
        if (data.parameters.contains(model_data_type::Parameters::Pvmin))
        {
          Pvmin_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Pvmin));
        }
        if (data.parameters.contains(model_data_type::Parameters::Pvmax))
        {
          Pvmax_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Pvmax));
        }
        if (data.parameters.contains(model_data_type::Parameters::T1))
        {
          T1_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::T1));
        }
        if (data.parameters.contains(model_data_type::Parameters::T2))
        {
          T2_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::T2));
        }
        if (data.parameters.contains(model_data_type::Parameters::T3))
        {
          T3_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::T3));
        }
        if (data.parameters.contains(model_data_type::Parameters::Dt))
        {
          Dt_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Dt));
        }
      }

      template <class ScalarT, typename IdxT>
      Tgov1<ScalarT, IdxT>::Tgov1(signal_type* pmech, signal_type* omega)
        : pmech_(pmech),
          omega_(omega),
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
        // Allocate local component data
        auto size = static_cast<size_t>(size_); // Avoid compiler warnings
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);

        // Set output signal after allocation
        // The signal is accessible to the generator
        if (pmech_)
        {
          pmech_->set(&y_[2]);
        }
        return 0;
      }

      /**
       * @brief Initialization of the Governor
       *
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::initialize()
      {
        ScalarT p0{0};

        // Initial mechanical = initial electric torque
        if (pmech_)
        {
          p0 = y_[2]; //<- generator needs to be initialized first
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
        ScalarT omega{0};
        if (omega_)
        {
          omega = omega_->read();
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
       * @brief The mechanical power output.
       * @warning This is not yet accessed by anything. The Genrou class will
       *          need to access this instead of a constant Pmech.
       * @return ScalarT - Mechanical output power value.
       */
      template <class ScalarT, typename IdxT>
      ScalarT& Tgov1<ScalarT, IdxT>::Pmech()
      {
        return y_[2];
      }

      // Available template instantiations
      template class Tgov1<double, long int>;
      template class Tgov1<double, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

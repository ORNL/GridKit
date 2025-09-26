#pragma once

/**
 * @file Tgov1Impl.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @brief Definition of a Turbine Governor Model (IEEET1).
 */

#include <cmath>
#include <iostream>

#include "Tgov1.hpp"
#include <Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#include <Model/PhasorDynamics/SignalNode/SignalNode.hpp>

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Constructs a Tgov1 governor model without setting its parameters
       *
       */
      template <class ScalarT, typename IdxT>
      Tgov1<ScalarT, IdxT>::Tgov1()
      {
        size_ = 3;
      }

      /**
       * @brief Constructs a Tgov1 governor model from its parameters
       *
       * @param pmech $P_m$ internal variable signal node
       * @param omega $\Delta_\omega$ external variable signal node
       */
      template <class ScalarT, typename IdxT>
      Tgov1<ScalarT, IdxT>::Tgov1(signal_type* pmech, signal_type* omega)
        : R_(0.05),
          Pvmin_(0),
          Pvmax_(1),
          T1_(0.5),
          T2_(2.5),
          T3_(7.5),
          Dt_(0)
      {
        signals_.template assignSignalNode<Tgov1InternalVariables::PM>(pmech);
        signals_.template attachSignalNode<Tgov1ExternalVariables::DELTAOMEGA>(omega);

        // 3 internal variables
        size_ = 3;
      }

      /**
       * @brief Constructs a Tgov1 governor model from its parameters
       *
       * @param data Data to initialize the model from.
       */
      template <class ScalarT, typename IdxT>
      Tgov1<ScalarT, IdxT>::Tgov1(const model_data_type& data)
      {
        initializeParameters(data);
        size_ = 3;
      }

      /**
       * @brief Helper function to extract and assign model parameters.
       *
       * Parses values from the model_data_type and assigns them to internal
       * parameters.
       *
       * @param data Structure containing model parameters.
       */
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

      /**
       * @brief Set the component ID
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /*!
       * @brief Allocate memory for model
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::allocate()
      {
        // Allocate local component data
        auto size = static_cast<size_t>(size_); // avoid compiler warnings
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);

        // Set output signal after allocation
        // The signal is accessible to the generator
        if (signals_.template isAssigned<Tgov1InternalVariables::PM>())
        {
          signals_.template getSignalNode<Tgov1InternalVariables::PM>()->set(&y_[2]);
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
        if (signals_.template isAssigned<Tgov1InternalVariables::PM>())
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
       * @brief Residuals of system equations evaluated locally
       *
       */
      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) int Tgov1<ScalarT, IdxT>::evaluateResidualLocally(ScalarT* y, ScalarT* yp, ScalarT* f)
      {
        // Read Internal Variables
        ScalarT ptx   = y[0]; // y0 - Ptx
        ScalarT pv    = y[1]; // y1 - Pv
        ScalarT pmech = y[2]; // y2 - Pmech

        // Read Internal Derivatives
        ScalarT ptx_dot = yp[0];
        ScalarT pv_dot  = yp[1];

        // The 'pre-limit' derivative of Pv
        ScalarT func     = (-pv + (pref_ - omega_) / R_) / T1_;
        ScalarT valv_ind = this->indicator(pv, func);

        // Internal Differential Equations
        f[0] = -ptx_dot + pv - (ptx + T2_ * pv) / T3_;
        f[1] = -pv_dot + valv_ind * func;

        // Internal Algebraic Equations
        f[2] = -pmech + (ptx + T2_ * pv) / T3_ - (Dt_ * omega_);

        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::evaluateResidual()
      {
        // Input Variables
        if (signals_.template isAttached<Tgov1ExternalVariables::DELTAOMEGA>())
        {
          omega_ = signals_.template readExternalVariable<Tgov1ExternalVariables::DELTAOMEGA>();
        }

        evaluateResidualLocally(y_.data(), yp_.data(), f_.data());

        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

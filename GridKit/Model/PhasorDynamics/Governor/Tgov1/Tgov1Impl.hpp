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

#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      using Log = ::GridKit::Utilities::Logger;

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
          T1_(HALF<RealT>),
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
          R_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::R));
        }

        if (data.parameters.contains(model_data_type::Parameters::Pvmin))
        {
          Pvmin_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Pvmin));
        }

        if (data.parameters.contains(model_data_type::Parameters::Pvmax))
        {
          Pvmax_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Pvmax));
        }

        if (data.parameters.contains(model_data_type::Parameters::T1))
        {
          T1_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::T1));
        }

        if (data.parameters.contains(model_data_type::Parameters::T2))
        {
          T2_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::T2));
        }

        if (data.parameters.contains(model_data_type::Parameters::T3))
        {
          T3_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::T3));
        }

        if (data.parameters.contains(model_data_type::Parameters::Dt))
        {
          Dt_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Dt));
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
        abs_tol_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        // Resize signal variable data
        ws_.resize(1);
        ws_indices_.resize(1);
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;

        // Default variable and residual index mapping to local index
        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        // Set output signals
        if (signals_.template isAssigned<Tgov1InternalVariables::PM>())
        {
          signals_.template getSignalNode<Tgov1InternalVariables::PM>()->set(&y_[2], &(this->getVariableIndex(2)));
        }

        return 0;
      }

      /**
       * @brief verify method checks that attached signals are also linked
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::verify() const
      {
        static constexpr auto DELTAOMEGA = Tgov1ExternalVariables::DELTAOMEGA;

        int ret = 0;

        if (signals_.template isAttached<DELTAOMEGA>())
        {
          if (!signals_.template isLinked<DELTAOMEGA>())
          {
            Log::error() << "Tgov1: deltaomega signal attached with no linked generator\n";
            ret += 1;
          }
        }

        return ret;
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
          p0 = y_[2]; ///<- generator needs to be initialized first
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
       * @brief Compute the absolute tolerance for each variable in the model
       * 
       * @param rel_tol The relative tolerance which can be used to pick the
       *        absolute tolerance.
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       * @return int 0 if successful, non-zero otherwise.
       * 
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
      {
        std::fill(abs_tol_.begin(), abs_tol_.end(), rel_tol);
        return 0;
      }

      /**
       * @brief Internal residuals
       *
       */
      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Tgov1<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT*                  y,
          ScalarT*                  yp,
          [[maybe_unused]] ScalarT* wb,
          ScalarT*                  ws,
          ScalarT*                  f)
      {
        // Read Internal Variables
        ScalarT ptx   = y[0]; // y0 - Ptx
        ScalarT pv    = y[1]; // y1 - Pv
        ScalarT pmech = y[2]; // y2 - Pmech

        // Read Internal Derivatives
        ScalarT ptx_dot = yp[0];
        ScalarT pv_dot  = yp[1];

        // Set signal variable aliases
        ScalarT omega = ws[0];

        // The 'pre-limit' derivative of Pv
        ScalarT func     = (-pv + (pref_ - omega) / R_) / T1_;
        ScalarT valv_ind = Math::indicator(Pvmin_, Pvmax_, pv, func);

        // Internal Differential Equations
        f[0] = -ptx_dot + pv - (ptx + T2_ * pv) / T3_;
        f[1] = -pv_dot + valv_ind * func;

        // Internal Algebraic Equations
        f[2] = -pmech + (ptx + T2_ * pv) / T3_ - (Dt_ * omega);

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
          ws_[0]         = signals_.template readExternalVariable<Tgov1ExternalVariables::DELTAOMEGA>();
          ws_indices_[0] = signals_.template readExternalVariableIndex<Tgov1ExternalVariables::DELTAOMEGA>();
        }

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

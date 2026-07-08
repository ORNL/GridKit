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
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
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
      template <typename scalar_type, typename index_type>
      Tgov1<scalar_type, index_type>::Tgov1()
        : Trate_(100.0)
      {
        size_ = 3;
        setDerivedParams();
      }

      /**
       * @brief Constructs a Tgov1 governor model from its parameters
       *
       * @param pmech $P_m$ internal variable signal node
       * @param omega $\Delta_\omega$ external variable signal node
       */
      template <typename scalar_type, typename index_type>
      Tgov1<scalar_type, index_type>::Tgov1(SignalNodeT* pmech, SignalNodeT* omega)
        : Trate_(100.0),
          R_(0.05),
          Pvmin_(0),
          Pvmax_(1),
          T1_(HALF<RealT>),
          T2_(2.5),
          T3_(7.5),
          Dt_(0)
      {
        ports_.out[Tgov1SignalOutputs::pmech].connect(pmech);
        ports_.in[Tgov1SignalInputs::speed].connect(omega);

        // 3 internal variables
        size_ = 3;
        setDerivedParams();
      }

      /**
       * @brief Constructs a Tgov1 governor model from its parameters
       *
       * @param data Data to initialize the model from.
       */
      template <typename scalar_type, typename index_type>
      Tgov1<scalar_type, index_type>::Tgov1(
          const ModelDataT& data, SignalNodeSetT& signal_nodes)
        : ports_(data, signal_nodes)
      {
        initializeParameters(data);

        size_ = 3;
        setDerivedParams();
      }

      /**
       * @brief Helper function to extract and assign model parameters.
       *
       * Parses values from the ModelDataT and assigns them to internal
       * parameters.
       *
       * @param data Structure containing model parameters.
       */
      template <typename scalar_type, typename index_type>
      void Tgov1<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;

        if (data.parameters.contains(Parameter::Trate))
        {
          Trate_ = std::get<RealT>(data.parameters.at(Parameter::Trate));
        }

        if (data.parameters.contains(Parameter::R))
        {
          R_ = std::get<RealT>(data.parameters.at(Parameter::R));
        }

        if (data.parameters.contains(Parameter::Pvmin))
        {
          Pvmin_ = std::get<RealT>(data.parameters.at(Parameter::Pvmin));
        }

        if (data.parameters.contains(Parameter::Pvmax))
        {
          Pvmax_ = std::get<RealT>(data.parameters.at(Parameter::Pvmax));
        }

        if (data.parameters.contains(Parameter::T1))
        {
          T1_ = std::get<RealT>(data.parameters.at(Parameter::T1));
        }

        if (data.parameters.contains(Parameter::T2))
        {
          T2_ = std::get<RealT>(data.parameters.at(Parameter::T2));
        }

        if (data.parameters.contains(Parameter::T3))
        {
          T3_ = std::get<RealT>(data.parameters.at(Parameter::T3));
        }

        if (data.parameters.contains(Parameter::Dt))
        {
          Dt_ = std::get<RealT>(data.parameters.at(Parameter::Dt));
        }
      }

      template <typename scalar_type, typename index_type>
      void Tgov1<scalar_type, index_type>::setDerivedParams()
      {
        va_component_base_ = Trate_ * static_cast<RealT>(1.0e6);
      }

      // System base -> component base when reading signals.
      template <typename scalar_type, typename index_type>
      scalar_type Tgov1<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * va_system_base_ / va_component_base_;
      }

      // Governor base -> system base for signals output.
      template <typename scalar_type, typename index_type>
      scalar_type Tgov1<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<scalar_type>(ONE<RealT>));
      }

      /**
       * @brief Set the component ID
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /*!
       * @brief Allocate memory for model
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::allocate()
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
        if (auto pmech_port = ports_.out[Tgov1SignalOutputs::pmech])
        {
          pmech_port.link(&y_[2], &(this->getVariableIndex(2)));
        }

        return 0;
      }

      /**
       * @brief verify method checks that attached signals are also linked
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::verify() const
      {
        int ret = 0;

        auto speed_port = ports_.in[Tgov1SignalInputs::speed];
        if (speed_port.connected() && !speed_port.linked())
        {
          Log::error() << "Tgov1: deltaomega signal attached with no linked generator\n";
          ret += 1;
        }

        return ret;
      }

      /**
       * @brief Initialization of the Governor
       *
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::initialize()
      {
        ScalarT p0{0};

        // Initial mechanical = initial electric torque
        if (ports_.out[Tgov1SignalOutputs::pmech])
        {
          // System base -> governor base for governor initialization.
          p0 = toComponentBase(y_[2]); ///<- generator needs to be initialized first
        }

        // Input Variables (Parameter for now)
        pref_ = R_ * p0;

        // Internal States
        y_[0] = (T3_ - T2_) * p0; // y0 - Ptx (Turbine Power )
        y_[1] = p0;               // y1 - Pv  (Valve Position)
        y_[2] = toSystemBase(p0); // y2 - Pm  (Mech Power, System Base)

        // D.V. Derivative
        yp_[0] = 0.0; // Ptx
        yp_[1] = 0.0; // Pv
        yp_[2] = 0.0; // Pm

        return 0;
      }

      /**
       * @brief Identify differential variables.
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::tagDifferentiable()
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
       * @tparam scalar_type Scalar data type
       * @tparam index_type Index data type
       * @return int 0 if successful, non-zero otherwise.
       *
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        std::fill(abs_tol_.begin(), abs_tol_.end(), rel_tol);
        return 0;
      }

      /**
       * @brief Internal residuals
       *
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int Tgov1<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
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

        // The 'pre-limit' target of Pv
        ScalarT func = -pv + (pref_ - omega) / R_;

        // Internal Differential Equations
        f[0] = -T3_ * ptx_dot - ptx + (T3_ - T2_) * pv;
        f[1] = -T1_ * pv_dot + Math::antiwindup(pv, func, Pvmin_, Pvmax_);

        // Internal Algebraic Equations
        // Convert pmech to component base from its system base value
        f[2] = -toComponentBase(pmech) + (ptx + T2_ * pv) / T3_ - (Dt_ * omega);

        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateResidual()
      {
        // Input Variables
        if (auto speed_port = ports_.in[Tgov1SignalInputs::speed])
        {
          ws_[0]         = speed_port.readSignal();
          ws_indices_[0] = speed_port.signalVariableIndex();
        }

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

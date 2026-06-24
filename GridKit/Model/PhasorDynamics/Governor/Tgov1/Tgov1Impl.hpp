#pragma once

/**
 * @file Tgov1Impl.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @brief Definition of a Turbine Governor Model (IEEET1).
 */

#include <iostream>

#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

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
      Tgov1<scalar_type, index_type>::Tgov1(SignalT* pmech, SignalT* omega)
        : Trate_(100.0),
          R_(0.05),
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
        setDerivedParams();
      }

      /**
       * @brief Constructs a Tgov1 governor model from its parameters
       *
       * @param data Data to initialize the model from.
       */
      template <typename scalar_type, typename index_type>
      Tgov1<scalar_type, index_type>::Tgov1(const ModelDataT& data)
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
        if (data.parameters.contains(ModelDataT::Parameters::Trate))
        {
          Trate_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Trate));
        }

        if (data.parameters.contains(ModelDataT::Parameters::R))
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
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        // Allocate local component data
        auto size = static_cast<size_t>(size_); // avoid compiler warnings

        tag_.resize(size);
        abs_tol_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        // Default variable and residual index mapping to local index
        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        // Resize signal variable data
        ws_.resize(1);
        ws_indices_.resize(1);
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;

        // Set output signals
        if (signals_.template isAssigned<Tgov1InternalVariables::PM>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<Tgov1InternalVariables::PM>()->set(&y[2], &(this->getVariableIndex(2)));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief verify method checks that attached signals are also linked
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::verify() const
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
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::initialize()
      {
        ScalarT p0{0};

        // Initial mechanical = initial electric torque
        auto* y  = y_.getData();
        auto* yp = yp_.getData();
        if (signals_.template isAssigned<Tgov1InternalVariables::PM>())
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
        yp[0] = 0.0; // Ptx
        yp[1] = 0.0; // Pv
        yp[2] = 0.0; // Pm

        y_.setDataUpdated();
        yp_.setDataUpdated();

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
      template <typename scalar_type, typename index_type>
      FORCE_INLINE int Tgov1<scalar_type, index_type>::evaluateInternalResidual(
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
        if (signals_.template isAttached<Tgov1ExternalVariables::DELTAOMEGA>())
        {
          ws_[0]         = signals_.template readExternalVariable<Tgov1ExternalVariables::DELTAOMEGA>();
          ws_indices_[0] = signals_.template readExternalVariableIndex<Tgov1ExternalVariables::DELTAOMEGA>();
        }

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);

        f_.setDataUpdated();

        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

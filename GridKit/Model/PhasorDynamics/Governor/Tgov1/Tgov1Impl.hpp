#pragma once

/**
 * @file Tgov1Impl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @brief Definition of the TGOV1 turbine-governor model.
 */

#include <algorithm>
#include <cmath>

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
        if (Pvmax_ < Pvmin_)
        {
          Log::warning() << "Tgov1: Pvmax is below Pvmin and the limits are swapped\n";
          std::swap(Pvmin_, Pvmax_);
        }
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

        variable_indices_.resize(size);
        residual_indices_.resize(size);

        // Default variable and residual index mapping to local index
        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        // Resize signal variable data
        const auto signal_size = static_cast<size_t>(Tgov1ExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        // Set output signals
        if (signals_.template isAssigned<Tgov1InternalVariables::PM>())
        {
          auto*      y  = y_.getData();
          const auto PM = static_cast<IdxT>(Tgov1InternalVariables::PM);
          signals_.template getSignalNode<Tgov1InternalVariables::PM>()->set(
              &y[static_cast<size_t>(PM)], &(this->getVariableIndex(PM)));
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
        static constexpr auto PREF       = Tgov1ExternalVariables::PREF;

        int ret = 0;

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Tgov1: " << message << '\n';
            ret += 1;
          }
        };

        check(std::isfinite(Trate_) && Trate_ > ZERO<RealT>,
              "Trate must be finite and positive");
        check(std::isfinite(va_system_base_) && va_system_base_ > ZERO<RealT>,
              "system power base must be finite and positive");
        check(std::isfinite(R_) && R_ != ZERO<RealT>, "R must be finite and nonzero");
        check(std::isfinite(Pvmin_) && std::isfinite(Pvmax_)
                  && std::isfinite(T1_) && std::isfinite(T2_)
                  && std::isfinite(T3_) && std::isfinite(Dt_),
              "parameters must be finite");
        check(signals_.template isAssigned<Tgov1InternalVariables::PM>(),
              "pmech output signal must be assigned");

        if (signals_.template isAttached<DELTAOMEGA>())
        {
          if (!signals_.template isLinked<DELTAOMEGA>())
          {
            Log::error() << "Tgov1: deltaomega signal attached with no linked generator\n";
            ret += 1;
          }
        }

        if (signals_.template isAttached<PREF>())
        {
          if (!signals_.template isLinked<PREF>())
          {
            Log::error() << "Tgov1: pref signal attached with no linked source\n";
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
        if (verify() != 0)
        {
          Log::error() << "Tgov1: cannot initialize with invalid configuration\n";
          return 1;
        }

        const auto PTX = static_cast<size_t>(Tgov1InternalVariables::PTX);
        const auto PV  = static_cast<size_t>(Tgov1InternalVariables::PV);
        const auto PM  = static_cast<size_t>(Tgov1InternalVariables::PM);

        auto* y = y_.getData();

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<Tgov1ExternalVariables::DELTAOMEGA>())
        {
          omega0 = signals_.template readExternalVariable<Tgov1ExternalVariables::DELTAOMEGA>();
        }

        const ScalarT pmech0 = y[PM];
        const ScalarT pm0    = toComponentBase(pmech0);
        const ScalarT pv0    = pm0 + Dt_ * omega0;
        const ScalarT pturb0 = pv0;
        const ScalarT pref0  = omega0 + R_ * pv0;

        if (!std::isfinite(static_cast<RealT>(pmech0))
            || !std::isfinite(static_cast<RealT>(omega0))
            || !std::isfinite(static_cast<RealT>(pv0))
            || !std::isfinite(static_cast<RealT>(pref0)))
        {
          Log::error() << "Tgov1: initial signals and states must be finite\n";
          return 1;
        }

        const RealT pv0_value = static_cast<RealT>(pv0);
        if (pv0_value < Pvmin_ || pv0_value > Pvmax_)
        {
          Log::warning() << "Tgov1: initial valve position is outside [Pvmin, Pvmax] "
                            "and the limits are widened\n";
        }
        Pvmin_ = std::min(Pvmin_, pv0_value);
        Pvmax_ = std::max(Pvmax_, pv0_value);

        y[PTX] = pturb0;
        y[PV]  = pv0;

        pref_set_ = pref0;
        if (signals_.template isAttached<Tgov1ExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<Tgov1ExternalVariables::PREF>(pref_set_);
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));

        return 0;
      }

      /**
       * @brief Identify differential variables.
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[static_cast<size_t>(Tgov1InternalVariables::PTX)] = T3_ != ZERO<RealT>;
        tag_[static_cast<size_t>(Tgov1InternalVariables::PV)]  = T1_ != ZERO<RealT>;
        tag_[static_cast<size_t>(Tgov1InternalVariables::PM)]  = false;

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
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
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
        const auto PTX = static_cast<size_t>(Tgov1InternalVariables::PTX);
        const auto PV  = static_cast<size_t>(Tgov1InternalVariables::PV);
        const auto PM  = static_cast<size_t>(Tgov1InternalVariables::PM);

        const auto DELTAOMEGA = static_cast<size_t>(Tgov1ExternalVariables::DELTAOMEGA);
        const auto PREF       = static_cast<size_t>(Tgov1ExternalVariables::PREF);

        const ScalarT pturb = y[PTX];
        const ScalarT pv    = y[PV];
        const ScalarT pmech = y[PM];

        const ScalarT pturb_dot = yp[PTX];
        const ScalarT pv_dot    = yp[PV];

        const ScalarT omega = ws[DELTAOMEGA];
        const ScalarT pref  = ws[PREF];

        const ScalarT valve_rate = -pv + (pref - omega) / R_;

        f[PTX] = -T3_ * pturb_dot - pturb + pv + T2_ * pv_dot;
        f[PV]  = -T1_ * pv_dot
                + Math::antiwindup(pv, valve_rate, Pvmin_, Pvmax_);
        f[PM] = -toComponentBase(pmech) + pturb - Dt_ * omega;

        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateResidual()
      {
        const auto DELTAOMEGA = static_cast<size_t>(Tgov1ExternalVariables::DELTAOMEGA);
        const auto PREF       = static_cast<size_t>(Tgov1ExternalVariables::PREF);

        ws_[DELTAOMEGA] = ScalarT{ZERO<RealT>};
        ws_[PREF]       = pref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<Tgov1ExternalVariables::DELTAOMEGA>())
        {
          ws_[DELTAOMEGA] = signals_.template readExternalVariable<Tgov1ExternalVariables::DELTAOMEGA>();
          ws_indices_[DELTAOMEGA] =
              signals_.template readExternalVariableIndex<Tgov1ExternalVariables::DELTAOMEGA>();
        }

        if (signals_.template isAttached<Tgov1ExternalVariables::PREF>())
        {
          ws_[PREF] = signals_.template readExternalVariable<Tgov1ExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<Tgov1ExternalVariables::PREF>();
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

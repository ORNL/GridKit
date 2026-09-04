#pragma once

/**
 * @file Tgov1Impl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the EMT TGOV1 turbine-governor model.
 */

#include <algorithm>
#include <limits>
#include <mutex>

#include <GridKit/Model/EMT/Component/Controller/TGOV1/Tgov1.hpp>
#include <GridKit/Model/EMT/Component/Controller/TGOV1/Tgov1Data.hpp>
#include <GridKit/Model/EMT/Signal/Signal.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
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
       * @param pmech $P_m$ internal variable signal
       * @param omega $\omega_r$ external variable signal
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
        signals_.template assignSignal<Tgov1InternalVariables::PM>(pmech);
        signals_.template attachSignal<Tgov1ExternalVariables::OMEGA>(omega);

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

      /**
       * @brief Static method to log time constant warnings
       *
       * @note Used in combination with static std:once_flag and std:call_once,
       *       to reduce the number of times the warning is printed.
       */
      template <typename scalar_type, typename index_type>
      void Tgov1<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "Tgov1: T1 and T3 below " << TIME_CONSTANT_MINIMUM
                       << " s are raised to that floor\n";
      }

      template <typename scalar_type, typename index_type>
      void Tgov1<scalar_type, index_type>::setDerivedParams()
      {
        if (T1_ < TIME_CONSTANT_MINIMUM || T3_ < TIME_CONSTANT_MINIMUM)
        {
          static std::once_flag time_constant_warning_flag_;
          std::call_once(time_constant_warning_flag_,
                         &logTimeConstantWarning);
        }

        T1_ = std::max(T1_, TIME_CONSTANT_MINIMUM);
        T3_ = std::max(T3_, TIME_CONSTANT_MINIMUM);
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

        // Resize coupling data
        this->allocateExternalVectors(static_cast<IdxT>(Tgov1ExternalVariables::MAXIMUM), 0);
        signals_.registerExternalVariableSignals(*this);

        // Set output signals
        if (signals_.template isAssigned<Tgov1InternalVariables::PM>())
        {
          auto*      y  = y_.getData();
          const auto PM = static_cast<IdxT>(Tgov1InternalVariables::PM);
          signals_.template getSignal<Tgov1InternalVariables::PM>()->set(
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
        static constexpr auto OMEGA = Tgov1ExternalVariables::OMEGA;
        static constexpr auto PREF  = Tgov1ExternalVariables::PREF;

        int ret = 0;

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Tgov1: " << message << '\n';
            ret += 1;
          }
        };

        check(Trate_ > ZERO<RealT>, "Trate must be positive");
        check(R_ != ZERO<RealT>, "R must be nonzero");
        check(Pvmin_ <= Pvmax_, "Pvmin must be less than or equal to Pvmax");
        check(signals_.template isAssigned<Tgov1InternalVariables::PM>(),
              "pmech output signal must be assigned");

        if (signals_.template isAttached<OMEGA>())
        {
          if (!signals_.template isLinked<OMEGA>())
          {
            Log::error() << "Tgov1: speed signal attached with no linked machine\n";
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

        ScalarT domega0{ZERO<RealT>};
        if (signals_.template isAttached<Tgov1ExternalVariables::OMEGA>())
        {
          domega0 = signals_.template readExternalVariable<Tgov1ExternalVariables::OMEGA>() - ONE<RealT>;
        }

        const ScalarT pm0    = y[PM];
        const ScalarT pv0    = pm0 + Dt_ * domega0;
        const ScalarT pturb0 = pv0;
        const ScalarT pref0  = domega0 + R_ * pv0;

        const RealT pv0_value       = static_cast<RealT>(pv0);
        const RealT limit_tolerance = static_cast<RealT>(4.0) * std::numeric_limits<RealT>::epsilon();
        if (pv0_value < Pvmin_ - limit_tolerance || pv0_value > Pvmax_ + limit_tolerance)
        {
          Log::error() << "Tgov1: initial valve position is outside [Pvmin, Pvmax]. "
                          "Check initial dispatch and valve limits\n";
          return 1;
        }

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
        tag_[0] = true;  // Ptx
        tag_[1] = true;  // Pv
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
          const ScalarT*                  y_ext,
          [[maybe_unused]] const ScalarT* yp_ext,
          ScalarT*                        f)
      {
        const auto PTX = static_cast<size_t>(Tgov1InternalVariables::PTX);
        const auto PV  = static_cast<size_t>(Tgov1InternalVariables::PV);
        const auto PM  = static_cast<size_t>(Tgov1InternalVariables::PM);

        const auto OMEGA = static_cast<size_t>(Tgov1ExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(Tgov1ExternalVariables::PREF);

        const ScalarT pturb = y[PTX];
        const ScalarT pv    = y[PV];
        const ScalarT pmech = y[PM];

        const ScalarT pturb_dot = yp[PTX];
        const ScalarT pv_dot    = yp[PV];

        const ScalarT domega = y_ext[OMEGA] - ONE<RealT>;
        const ScalarT pref   = y_ext[PREF];

        f[PTX] = -pturb_dot - (pturb - pv - T2_ * pv_dot) / T3_;
        f[PV]  = -pv_dot + Math::antiwindup(pv, -pv + (pref - domega) / R_, Pvmin_, Pvmax_) / T1_;
        f[PM]  = -pmech + pturb - Dt_ * domega;

        return 0;
      }

      /**
       * @brief Gather external variables and index maps.
       *
       * The synchronous-speed default backs the speed input and the latched
       * reference backs the governor reference when no producer is attached;
       * the base gather then refreshes every attached slot through its
       * signal.
       */
      template <typename scalar_type, typename index_type>
      void Tgov1<scalar_type, index_type>::gatherExternalVariables()
      {
        y_ext_[static_cast<size_t>(Tgov1ExternalVariables::OMEGA)] = ScalarT{ONE<RealT>};
        y_ext_[static_cast<size_t>(Tgov1ExternalVariables::PREF)]  = pref_set_;

        Component<scalar_type, index_type>::gatherExternalVariables();
      }

      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateInternalResidual()
      {
        gatherExternalVariables();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f);

        f_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateResidual()
      {
        evaluateInternalResidual();
        return this->evaluateExternalResidual();
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

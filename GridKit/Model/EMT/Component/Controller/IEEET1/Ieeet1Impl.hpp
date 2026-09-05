/**
 * @file Ieeet1Impl.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Definition of a IEEET1 Exciter.
 *
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <mutex>

#pragma once

#include <limits>

#include <GridKit/Model/EMT/Component/Controller/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/EMT/Component/Controller/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/EMT/Signal/Signal.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief  Constructor for IEEET1 Exciter
       */
      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::Ieeet1()
        : Ieeet1(ModelDataT{})
      {
      }

      /**
       * @brief  Constructor for IEEET1 Exciter
       *
       * @param data  Data object to store parameters
       */
      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::Ieeet1(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        // Parse data struct into model
        this->initModelParams(data);

        initializeMonitor();

        // 9 Internal Variables
        size_ = 9;
      }

      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::~Ieeet1()
      {
      }

      /**
       * @brief Set the component ID
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate memory for model
       *
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
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

        this->allocateExternalVectors(static_cast<IdxT>(Ieeet1ExternalVariables::MAXIMUM), 0);
        signals_.registerExternalVariableSignals(*this);
        uel_on_ = signals_.template isAttached<Ieeet1ExternalVariables::VUEL>() ? ONE<RealT> : ZERO<RealT>;
        oel_on_ = signals_.template isAttached<Ieeet1ExternalVariables::VOEL>() ? ONE<RealT> : ZERO<RealT>;

        // Set output signals
        if (signals_.template isAssigned<Ieeet1InternalVariables::EFD>())
        {
          auto* y = y_.getData();
          signals_.template getSignal<Ieeet1InternalVariables::EFD>()->set(&y[7], &(this->getVariableIndex(7)));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Verify parameter values and attached signal links
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::verify() const
      {
        int ret = 0;

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Ieeet1: " << message << '\n';
            ret += 1;
          }
        };

        check(V_ > ZERO<RealT>, "V must be positive");
        for (const RealT value : {V_, Tr_, Ka_, Ta_, Ke_, Te_, Kf_, Tf_, Vrmin_, Vrmax_, E1_, E2_, Se1_, Se2_, Ispdlim_, SA_, SB_})
        {
          check(std::isfinite(value), "parameters and derived saturation coefficients must be finite");
        }
        check(signals_.template isAssigned<Ieeet1InternalVariables::EFD>(), "efd output signal must be assigned");
        check(signals_.template isAttached<Ieeet1ExternalVariables::VA>()
                  && signals_.template isAttached<Ieeet1ExternalVariables::VB>()
                  && signals_.template isAttached<Ieeet1ExternalVariables::VC>(),
              "bus voltage port must be attached");
        check(Ka_ > ZERO<RealT>, "Ka must be positive");
        check(Vrmin_ <= Vrmax_, "Vrmin must be less than or equal to Vrmax");
        check(Ispdlim_ == ZERO<RealT> || Ispdlim_ == ONE<RealT>,
              "Ispdlim must be 0 or 1");

        const bool saturation_disabled =
            Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>;

        if (!saturation_disabled)
        {
          check(E1_ > ZERO<RealT>, "E1 must be positive when saturation is enabled");
          check(E2_ > ZERO<RealT>, "E2 must be positive when saturation is enabled");
          check(Se1_ >= ZERO<RealT>, "Se1 must be non-negative when saturation is enabled");
          check(Se2_ >= ZERO<RealT>, "Se2 must be non-negative when saturation is enabled");

          const bool sat_ordered = (E2_ > E1_ && Se2_ > Se1_) || (E2_ < E1_ && Se2_ < Se1_);
          check(sat_ordered, "E1/E2 and Se1/Se2 must be ordered consistently");
        }

        auto check_attached_signal =
            [&]<Ieeet1ExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "Ieeet1: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_attached_signal.template operator()<Ieeet1ExternalVariables::VA>("bus phase a");
        check_attached_signal.template operator()<Ieeet1ExternalVariables::VB>("bus phase b");
        check_attached_signal.template operator()<Ieeet1ExternalVariables::VC>("bus phase c");
        check_attached_signal.template operator()<Ieeet1ExternalVariables::OMEGA>("speed");
        check_attached_signal.template operator()<Ieeet1ExternalVariables::VREF>("vref");
        check_attached_signal.template operator()<Ieeet1ExternalVariables::VS>("vs");
        check_attached_signal.template operator()<Ieeet1ExternalVariables::VUEL>("vuel");
        check_attached_signal.template operator()<Ieeet1ExternalVariables::VOEL>("voel");

        return ret;
      }

      /** Initialize after the machine has seeded the field-voltage output. */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::initialize()
      {
        if (!allocated_ || verify() != 0)
        {
          Log::error() << "Ieeet1: cannot initialize with invalid configuration\n";
          return 1;
        }
        gatherExternalVariables();
        const auto*   external     = y_ext_.data();
        const ScalarT omega        = external[3];
        const ScalarT vs           = external[5];
        const ScalarT vuel         = external[6];
        const ScalarT voel         = external[7];
        const ScalarT speed_factor = ONE<RealT> + (omega - ONE<RealT>) *Ispdlim_;
        auto*         y            = y_.getData();
        const ScalarT efd0         = y[7];
        if (!(static_cast<RealT>(speed_factor) > ZERO<RealT>)
            || !std::isfinite(static_cast<RealT>(speed_factor)))
        {
          Log::error() << "Ieeet1: initial field-voltage speed multiplier must be finite and positive\n";
          return 1;
        }
        const ScalarT efdp   = efd0 / speed_factor;
        const ScalarT ksat   = SB_ * Math::qramp(efdp - SA_);
        RealT         ke_eff = Ke_;
        if (Ke_ == ZERO<RealT>)
        {
          if (static_cast<RealT>(efdp) == ZERO<RealT>)
          {
            Log::error() << "Ieeet1: automatic Ke requires nonzero initial field voltage\n";
            return 1;
          }
          ke_eff = (Vrmax_ / static_cast<RealT>(10) - static_cast<RealT>(ksat)) / static_cast<RealT>(efdp);
        }
        const ScalarT ec   = voltageMagnitude(external);
        const ScalarT vr   = ke_eff * efdp + ksat;
        const ScalarT vtr  = vr / Ka_;
        const ScalarT vfx  = (Kf_ / Tf_) * efdp;
        const ScalarT vref = ec + vtr - vs - vuel - voel;
        for (const ScalarT& value : {efd0, efdp, ksat, ec, vr, vtr, vfx, vref})
        {
          if (!std::isfinite(static_cast<RealT>(value)) || !std::isfinite(ke_eff))
          {
            Log::error() << "Ieeet1: initial states and effective Ke must be finite\n";
            return 1;
          }
        }
        const RealT tolerance = static_cast<RealT>(4) * std::numeric_limits<RealT>::epsilon()
                                * std::max({ONE<RealT>, std::abs(Vrmin_), std::abs(Vrmax_)});
        if (static_cast<RealT>(vr) < Vrmin_ - tolerance || static_cast<RealT>(vr) > Vrmax_ + tolerance)
        {
          Log::error() << "Ieeet1: initial regulator voltage is outside [Vrmin, Vrmax]\n";
          return 1;
        }

        Ke_eff_   = ke_eff;
        y[0]      = ec;
        y[1]      = vr;
        y[2]      = efdp;
        y[3]      = vfx;
        y[4]      = vtr;
        y[5]      = ZERO<RealT>;
        y[6]      = ksat;
        y[8]      = ksat;
        vref_set_ = vref;
        if (signals_.template isAttached<Ieeet1ExternalVariables::VREF>())
        {
          signals_.template writeExternalVariable<Ieeet1ExternalVariables::VREF>(vref);
        }
        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief  Identify differential variables.
       *
       * @return int 0
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[0] = true;  // y0 - vts  - Sensed term volt
        tag_[1] = true;  // y1 - vr   - Voltage reg
        tag_[2] = true;  // y2 - efdp - Efd pre mult
        tag_[3] = true;  // y3 - vfx  - Exciter feedback
        tag_[4] = false; // y4 - vtr  - Term Volt Err
        tag_[5] = false; // y5 - vf   - Feedback volt
        tag_[6] = false; // y6 - ve   - Excit. Cntrl Volt
        tag_[7] = false; // y7 - efd  - Efd
        tag_[8] = false; // y8 - ksat - Saturation

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
      int Ieeet1<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Internal Residual
       *
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int Ieeet1<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          const ScalarT*                  external,
          [[maybe_unused]] const ScalarT* external_dot,
          ScalarT*                        f)
      {
        const auto OMEGA = static_cast<size_t>(Ieeet1ExternalVariables::OMEGA);
        const auto VREF  = static_cast<size_t>(Ieeet1ExternalVariables::VREF);
        const auto VS    = static_cast<size_t>(Ieeet1ExternalVariables::VS);
        const auto VUEL  = static_cast<size_t>(Ieeet1ExternalVariables::VUEL);
        const auto VOEL  = static_cast<size_t>(Ieeet1ExternalVariables::VOEL);

        // Read bus voltage components
        ScalarT Ec = voltageMagnitude(external);

        // Read Internal Variables
        ScalarT vts  = y[0]; // y0 - Sensed term volt
        ScalarT vr   = y[1]; // y1 - Voltage reg
        ScalarT efdp = y[2]; // y2 - Efd pre mult
        ScalarT vfx  = y[3]; // y3 - Exciter feedback
        ScalarT vtr  = y[4]; // y4 - Term Volt Err
        ScalarT vf   = y[5]; // y5 - Feedback volt
        ScalarT ve   = y[6]; // y6 - Excit. Cntrl Volt
        ScalarT efd  = y[7]; // y7 - Efd
        ScalarT ksat = y[8]; // y8 - Saturation

        // Read Internal Derivatives
        ScalarT vts_dot  = yp[0];
        ScalarT vr_dot   = yp[1];
        ScalarT efdp_dot = yp[2];
        ScalarT vfx_dot  = yp[3];

        // Set signal variable aliases
        ScalarT omega = external[OMEGA];
        ScalarT vref  = external[VREF];
        ScalarT vs    = external[VS];
        ScalarT vuel  = external[VUEL];
        ScalarT voel  = external[VOEL];

        // The 'pre-limit' derivative of Vr.
        ScalarT func = (-vr + Ka_ * vtr) / Ta_;

        // Internal Differential Equations
        f[0] = -vts_dot + (Ec - vts) / Tr_;
        f[1] = -vr_dot + Math::antiwindup(vr, func, Vrmin_, Vrmax_);
        f[2] = -efdp_dot + (vr - ve - Ke_eff_ * efdp) / Te_;
        f[3] = -vfx_dot + vf / Tf_;

        // Internal Algebraic Equations
        f[4] = -vts + vref + vs + uel_on_ * vuel + oel_on_ * voel - vtr - vf;
        f[5] = -Tf_ * (vf + vfx) + Kf_ * efdp;
        f[6] = -ve + ksat;
        f[7] = -efd + efdp + (omega - ONE<RealT>) *efdp * Ispdlim_;
        f[8] = -ksat + SB_ * Math::qramp(efdp - SA_);

        return 0;
      }

      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::gatherExternalVariables()
      {
        y_ext_[3] = omega_set_;
        y_ext_[4] = vref_set_;
        y_ext_[5] = vs_set_;
        y_ext_[6] = vuel_set_;
        y_ext_[7] = voel_set_;
        Component<scalar_type, index_type>::gatherExternalVariables();
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateInternalResidual()
      {
        gatherExternalVariables();
        evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      /**
       * @brief Initialization Exciter Parameters from data structure
       */
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;

        if (data.parameters.contains(Parameter::V))
        {
          V_ = std::get<RealT>(data.parameters.at(Parameter::V));
        }
        if (data.parameters.contains(Parameter::Tr))
        {
          Tr_ = std::get<RealT>(data.parameters.at(Parameter::Tr));
        }
        if (data.parameters.contains(Parameter::Ka))
        {
          Ka_ = std::get<RealT>(data.parameters.at(Parameter::Ka));
        }
        if (data.parameters.contains(Parameter::Ta))
        {
          Ta_ = std::get<RealT>(data.parameters.at(Parameter::Ta));
        }
        if (data.parameters.contains(Parameter::Ke))
        {
          Ke_ = std::get<RealT>(data.parameters.at(Parameter::Ke));
        }
        if (data.parameters.contains(Parameter::Te))
        {
          Te_ = std::get<RealT>(data.parameters.at(Parameter::Te));
        }
        if (data.parameters.contains(Parameter::Kf))
        {
          Kf_ = std::get<RealT>(data.parameters.at(Parameter::Kf));
        }
        if (data.parameters.contains(Parameter::Tf))
        {
          Tf_ = std::get<RealT>(data.parameters.at(Parameter::Tf));
        }
        if (data.parameters.contains(Parameter::Vrmin))
        {
          Vrmin_ = std::get<RealT>(data.parameters.at(Parameter::Vrmin));
        }
        if (data.parameters.contains(Parameter::Vrmax))
        {
          Vrmax_ = std::get<RealT>(data.parameters.at(Parameter::Vrmax));
        }
        if (data.parameters.contains(Parameter::E1))
        {
          E1_ = std::get<RealT>(data.parameters.at(Parameter::E1));
        }
        if (data.parameters.contains(Parameter::E2))
        {
          E2_ = std::get<RealT>(data.parameters.at(Parameter::E2));
        }
        if (data.parameters.contains(Parameter::Se1))
        {
          Se1_ = std::get<RealT>(data.parameters.at(Parameter::Se1));
        }
        if (data.parameters.contains(Parameter::Se2))
        {
          Se2_ = std::get<RealT>(data.parameters.at(Parameter::Se2));
        }
        if (data.parameters.contains(Parameter::Ispdlim))
        {
          Ispdlim_ = std::get<RealT>(data.parameters.at(Parameter::Ispdlim));
        }

        Ke_eff_ = Ke_;
        setDerivedParameters();
      }

      /**
       * @brief Static method to log time constant warnings
       *
       * @note Used in combination with static std:once_flag and std:call_once,
       *       to reduce the number of times the warning is printed.
       */
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "Ieeet1: Tr, Ta, Te, and Tf below "
                       << TIME_CONSTANT_MINIMUM
                       << " s are raised to that floor\n";
      }

      /**
       * @brief Resolve the parameter-derived constants
       */
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::setDerivedParameters()
      {
        for (const RealT value : {Tr_, Ta_, Te_, Tf_})
        {
          if (!std::isfinite(value))
            return;
        }
        if (Tr_ < TIME_CONSTANT_MINIMUM || Ta_ < TIME_CONSTANT_MINIMUM
            || Te_ < TIME_CONSTANT_MINIMUM || Tf_ < TIME_CONSTANT_MINIMUM)
        {
          static std::once_flag time_constant_warning_flag_;
          std::call_once(time_constant_warning_flag_,
                         &logTimeConstantWarning);
        }

        Tr_ = std::max(Tr_, TIME_CONSTANT_MINIMUM);
        Ta_ = std::max(Ta_, TIME_CONSTANT_MINIMUM);
        Te_ = std::max(Te_, TIME_CONSTANT_MINIMUM);
        Tf_ = std::max(Tf_, TIME_CONSTANT_MINIMUM);

        SA_ = ZERO<RealT>;
        SB_ = ZERO<RealT>;

        const bool saturation_disabled =
            Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>;

        if (saturation_disabled)
        {
          return;
        }

        const bool sat_ordered = (E2_ > E1_ && Se2_ > Se1_) || (E2_ < E1_ && Se2_ < Se1_);
        if (E1_ <= ZERO<RealT> || E2_ <= ZERO<RealT>
            || Se1_ < ZERO<RealT> || Se2_ < ZERO<RealT>
            || !sat_ordered)
        {
          return;
        }

        if (Se1_ == ZERO<RealT>)
        {
          const RealT dE = E2_ - E1_;
          SA_            = E1_;
          SB_            = Se2_ * E2_ / (dE * dE);
          return;
        }

        if (Se2_ == ZERO<RealT>)
        {
          const RealT dE = E1_ - E2_;
          SA_            = E2_;
          SB_            = Se1_ * E1_ / (dE * dE);
          return;
        }

        const RealT C = std::sqrt(Se2_ * E2_ / (Se1_ * E1_));

        // Solution 1 (Aligned with PW)
        SA_ = (C * E1_ - E2_) / (C - ONE<RealT>);
        SB_ = Se1_ * E1_ / ((E1_ - SA_) * (E1_ - SA_));
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Ieeet1<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = ModelDataT::MonitorableVariables;
        monitor_->set(Variable::efd, [this]
                      { return y_.getData()[7]; });
        monitor_->set(Variable::ksat, [this]
                      { return y_.getData()[8]; });
        monitor_->set(Variable::vts, [this]
                      { return y_.getData()[0]; });
        monitor_->set(Variable::vr, [this]
                      { return y_.getData()[1]; });
        monitor_->set(Variable::vref, [this]
                      {
                        if (signals_.template isAttached<Ieeet1ExternalVariables::VREF>())
                          return signals_.template readExternalVariable<Ieeet1ExternalVariables::VREF>();
                        return vref_set_; });
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

/**
 * @file Esdc1aImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the ESDC1A exciter model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1a.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Logger used for ESDC1A diagnostics.
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct an ESDC1A exciter without parameters
       *
       * The model is sized but left unconfigured. Every parameter keeps its
       * documented default, and no monitor is created. verify() rejects the
       * model until an `efd` output node is assigned.
       *
       * @param[in] bus Terminal bus the exciter measures.
       */
      template <typename scalar_type, typename index_type>
      Esdc1a<scalar_type, index_type>::Esdc1a(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(Esdc1aInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      /**
       * @brief Construct an ESDC1A exciter from model data
       *
       * @param[in] bus Terminal bus the exciter measures.
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      Esdc1a<scalar_type, index_type>::Esdc1a(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(Esdc1aInternalVariables::MAXIMUM);
      }

      /**
       * @brief Destroy the ESDC1A exciter.
       */
      template <typename scalar_type, typename index_type>
      Esdc1a<scalar_type, index_type>::~Esdc1a()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the model vectors and wire the field-voltage output
       *
       * Sizes the state, residual, bus-interface, and signal-interface
       * buffers, seeds the identity index maps, and points an assigned `efd`
       * node at the internal field-voltage state. That node aliases ESDC1A
       * storage from here on, which is how initialize() reads the seed a
       * machine model wrote. Repeated calls reuse the allocated vectors.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::allocate()
      {
        const auto EFD = static_cast<size_t>(Esdc1aInternalVariables::EFD);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(2, ScalarT{0});

        const auto signal_size = static_cast<size_t>(Esdc1aExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        auto* y = y_.getData();

        if (signals_.template isAssigned<Esdc1aInternalVariables::EFD>())
        {
          signals_.template getSignalNode<Esdc1aInternalVariables::EFD>()->set(
              &y[EFD],
              &(this->getVariableIndex(static_cast<IdxT>(EFD))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the ESDC1A configuration
       *
       * Checks parameter-loading errors, static parameter relationships,
       * terminal-bus association, the required field-voltage output, and
       * attached external signals. Seed feasibility is operating-point
       * dependent and is checked by initialize().
       *
       * @return Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Esdc1a: " << message << '\n';
            ret += 1;
          }
        };

        if (bus_ == nullptr)
        {
          Log::error() << "Esdc1a: bus pointer is null\n";
          ret += 1;
        }

        check(Ka_ > ZERO<RealT>, "Ka must be positive");
        check(Tc_ >= ZERO<RealT>, "Tc must be non-negative");
        check(Vrmin_ <= Vrmax_, "Vrmin must be less than or equal to Vrmax");
        check(UEL_ >= static_cast<IdxT>(0) && UEL_ <= static_cast<IdxT>(3),
              "UEL must be 0, 1, 2, or 3");

        if (!(Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>) )
        {
          check(E1_ > ZERO<RealT>, "E1 must be positive when saturation is enabled");
          check(E2_ > ZERO<RealT>, "E2 must be positive when saturation is enabled");
          check(Se1_ > ZERO<RealT>, "Se1 must be positive when saturation is enabled");
          check(Se2_ > ZERO<RealT>, "Se2 must be positive when saturation is enabled");

          const bool saturation_points_are_ordered =
              (E2_ > E1_ && Se2_ > Se1_)
              || (E2_ < E1_ && Se2_ < Se1_);
          check(saturation_points_are_ordered,
                "E1/E2 and Se1/Se2 must be ordered consistently");
        }

        if (!signals_.template isAssigned<Esdc1aInternalVariables::EFD>())
        {
          Log::error() << "Esdc1a: required efd output signal is not assigned\n";
          ret += 1;
        }

        if (Spdmlt_ && !signals_.template isAttached<Esdc1aExternalVariables::OMEGA>())
        {
          Log::error() << "Esdc1a: speed signal is required when Spdmlt is enabled\n";
          ret += 1;
        }

        // An attached port must resolve to writable signal storage. The
        // enumerator is a template argument, so each port names itself once.
        auto check_attached_signal =
            [&]<Esdc1aExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "Esdc1a: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_attached_signal.template operator()<Esdc1aExternalVariables::OMEGA>("speed");
        check_attached_signal.template operator()<Esdc1aExternalVariables::VREF>("vref");
        check_attached_signal.template operator()<Esdc1aExternalVariables::VS>("vs");
        check_attached_signal.template operator()<Esdc1aExternalVariables::VUEL>("vuel");

        return ret;
      }

      /**
       * @brief Initialize ESDC1A from the field-voltage output
       *
       * Resolves the steady internal state and voltage reference while
       * preserving the seeded `efd`, latches attached Known inputs, and
       * publishes the reference to an attached `vref` signal.
       *
       * @return Zero on success; nonzero when the configuration or operating point is rejected.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::initialize()
      {
        const auto EFDP = static_cast<size_t>(Esdc1aInternalVariables::EFDP);
        const auto VC   = static_cast<size_t>(Esdc1aInternalVariables::VC);
        const auto VR   = static_cast<size_t>(Esdc1aInternalVariables::VR);
        const auto VF   = static_cast<size_t>(Esdc1aInternalVariables::VF);
        const auto XLL  = static_cast<size_t>(Esdc1aInternalVariables::XLL);
        const auto EV   = static_cast<size_t>(Esdc1aInternalVariables::EV);
        const auto VLL  = static_cast<size_t>(Esdc1aInternalVariables::VLL);
        const auto VHV  = static_cast<size_t>(Esdc1aInternalVariables::VHV);
        const auto SE   = static_cast<size_t>(Esdc1aInternalVariables::SE);
        const auto VFE  = static_cast<size_t>(Esdc1aInternalVariables::VFE);
        const auto EFD  = static_cast<size_t>(Esdc1aInternalVariables::EFD);

        bool ret = verify() == 0;
        if (!ret)
        {
          Log::error() << "Esdc1a: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y = y_.getData();

        // The assigned efd node aliases this entry after allocate(). Its
        // seeded value remains untouched throughout initialization.
        const ScalarT efd0 = y[EFD];

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<Esdc1aExternalVariables::OMEGA>())
        {
          omega0 = signals_.template readExternalVariable<Esdc1aExternalVariables::OMEGA>();
        }

        ScalarT vs0{ZERO<RealT>};
        if (signals_.template isAttached<Esdc1aExternalVariables::VS>())
        {
          vs0 = signals_.template readExternalVariable<Esdc1aExternalVariables::VS>();
        }

        ScalarT vuel0{ZERO<RealT>};
        if (signals_.template isAttached<Esdc1aExternalVariables::VUEL>())
        {
          vuel0 = signals_.template readExternalVariable<Esdc1aExternalVariables::VUEL>();
        }

        const ScalarT vc0 = std::sqrt(Vr() * Vr() + Vi() * Vi());

        ret = std::isfinite(static_cast<RealT>(efd0))
              && std::isfinite(static_cast<RealT>(vc0));
        if (!ret)
        {
          Log::error() << "Esdc1a: initial bus voltage and field-voltage seed must be finite\n";
          return 1;
        }

        const ScalarT speed_multiplier = ONE<RealT> + spd_on_ * omega0;

        ret = std::isfinite(static_cast<RealT>(speed_multiplier))
              && speed_multiplier > ZERO<RealT>;
        if (!ret)
        {
          Log::error() << "Esdc1a: speed multiplier must be finite and positive at initialization\n";
          return 1;
        }

        const ScalarT efdp0 = efd0 / speed_multiplier;

        ret = !exclim_ || efdp0 >= ZERO<RealT>;
        if (!ret)
        {
          Log::error() << "Esdc1a: initial Efd' is below its enabled zero limit\n";
          return 1;
        }

        const ScalarT se0  = SB_ * Math::qramp(efdp0 - SA_);
        const ScalarT vfe0 = (Ke_ + se0) * efdp0;
        const ScalarT vr0  = vfe0;
        const ScalarT vhv0 = vr0 / Ka_;

        ret = vr0 >= Vrmin_ && vr0 <= Vrmax_;
        if (!ret)
        {
          Log::error() << "Esdc1a: initialized VR is outside limits\n";
          return 1;
        }

        // An inactive high-value gate is seeded with the gate input, so the
        // residual reproduces VHV through the same smooth maximum.
        ScalarT vll0 = vhv0;
        if (uel_on_ == ZERO<RealT>)
        {
          const RealT gate_margin0 = static_cast<RealT>(vhv0 - vuel0);

          ret = gate_margin0 > ZERO<RealT>;
          if (!ret)
          {
            Log::error() << "Esdc1a: smooth high-value gate is active at initialization\n";
            return 1;
          }
          vll0 = vuel0 + inverseRamp(gate_margin0);
        }

        const ScalarT vf0   = ScalarT{ZERO<RealT>};
        const ScalarT ev0   = vll0;
        const ScalarT xll0  = ev0;
        const ScalarT vref0 = ev0 + vc0 + vf0 - vs0 - uel_on_ * vuel0;

        y[EFDP] = efdp0;
        y[VC]   = vc0;
        y[VR]   = vr0;
        y[VF]   = vf0;
        y[XLL]  = xll0;
        y[EV]   = ev0;
        y[VLL]  = vll0;
        y[VHV]  = vhv0;
        y[SE]   = se0;
        y[VFE]  = vfe0;
        y[EFD]  = efd0;

        omega_set_ = omega0;
        vref_set_  = vref0;
        vs_set_    = vs0;
        vuel_set_  = vuel0;

        if (signals_.template isAttached<Esdc1aExternalVariables::VREF>())
        {
          signals_.template writeExternalVariable<Esdc1aExternalVariables::VREF>(vref_set_);
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief Identify the differential variables
       *
       * The field-voltage state, the voltage transducer, the regulator, the
       * stabilizing feedback, and the lead-lag state carry derivatives;
       * every other internal variable is algebraic.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::tagDifferentiable()
      {
        const auto EFDP = static_cast<size_t>(Esdc1aInternalVariables::EFDP);
        const auto VC   = static_cast<size_t>(Esdc1aInternalVariables::VC);
        const auto VR   = static_cast<size_t>(Esdc1aInternalVariables::VR);
        const auto VF   = static_cast<size_t>(Esdc1aInternalVariables::VF);
        const auto XLL  = static_cast<size_t>(Esdc1aInternalVariables::XLL);

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[EFDP] = true;
        tag_[VC]   = true;
        tag_[VR]   = true;
        tag_[VF]   = true;
        tag_[XLL]  = true;
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * Every internal variable receives @p rel_tol as its absolute
       * tolerance.
       *
       * @param[in] rel_tol Solver relative tolerance.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       * Refreshes the bus and signal interface buffers and evaluates the
       * internal residual. ESDC1A injects no current, so there is no bus
       * residual. An unattached input port falls back to the value latched
       * by initialize().
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::evaluateResidual()
      {
        const auto OMEGA = static_cast<size_t>(Esdc1aExternalVariables::OMEGA);
        const auto VREF  = static_cast<size_t>(Esdc1aExternalVariables::VREF);
        const auto VS    = static_cast<size_t>(Esdc1aExternalVariables::VS);
        const auto VUEL  = static_cast<size_t>(Esdc1aExternalVariables::VUEL);

        ws_[OMEGA] = omega_set_;
        ws_[VREF]  = vref_set_;
        ws_[VS]    = vs_set_;
        ws_[VUEL]  = vuel_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<Esdc1aExternalVariables::OMEGA>())
        {
          ws_[OMEGA] = signals_.template readExternalVariable<Esdc1aExternalVariables::OMEGA>();
          ws_indices_[OMEGA] =
              signals_.template readExternalVariableIndex<Esdc1aExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<Esdc1aExternalVariables::VREF>())
        {
          ws_[VREF] = signals_.template readExternalVariable<Esdc1aExternalVariables::VREF>();
          ws_indices_[VREF] =
              signals_.template readExternalVariableIndex<Esdc1aExternalVariables::VREF>();
        }
        if (signals_.template isAttached<Esdc1aExternalVariables::VS>())
        {
          ws_[VS] = signals_.template readExternalVariable<Esdc1aExternalVariables::VS>();
          ws_indices_[VS] =
              signals_.template readExternalVariableIndex<Esdc1aExternalVariables::VS>();
        }
        if (signals_.template isAttached<Esdc1aExternalVariables::VUEL>())
        {
          ws_[VUEL] = signals_.template readExternalVariable<Esdc1aExternalVariables::VUEL>();
          ws_indices_[VUEL] =
              signals_.template readExternalVariableIndex<Esdc1aExternalVariables::VUEL>();
        }

        wb_[0] = Vr();
        wb_[1] = Vi();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);
        f_.setDataUpdated();
        return 0;
      }

      /**
       * @brief Access the monitor
       *
       * @return Monitor for this model, or nullptr when the model was
       *         constructed without data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Esdc1a<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Evaluate the ESDC1A internal residual.
       *
       * Evaluates the five exciter states and the six algebraic equations
       * documented in the model README. The body is kept free of branches
       * and loops so sparse automatic differentiation resolves a fixed
       * structure; the three selector decisions enter as multiplicative
       * masks set by setDerivedParameters().
       *
       * @param[in] y Internal variables in Esdc1aInternalVariables order.
       * @param[in] yp Internal derivatives in the same enum order.
       * @param[in] wb Terminal-bus \f$(V_{\mathrm{r}},V_{\mathrm{i}})\f$
       *               voltage components.
       * @param[in] ws Signal values in Esdc1aExternalVariables order.
       * @param[out] f Residuals in Esdc1aInternalVariables order.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Esdc1a<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto EFDP = static_cast<size_t>(Esdc1aInternalVariables::EFDP);
        const auto VC   = static_cast<size_t>(Esdc1aInternalVariables::VC);
        const auto VR   = static_cast<size_t>(Esdc1aInternalVariables::VR);
        const auto VF   = static_cast<size_t>(Esdc1aInternalVariables::VF);
        const auto XLL  = static_cast<size_t>(Esdc1aInternalVariables::XLL);
        const auto EV   = static_cast<size_t>(Esdc1aInternalVariables::EV);
        const auto VLL  = static_cast<size_t>(Esdc1aInternalVariables::VLL);
        const auto VHV  = static_cast<size_t>(Esdc1aInternalVariables::VHV);
        const auto SE   = static_cast<size_t>(Esdc1aInternalVariables::SE);
        const auto VFE  = static_cast<size_t>(Esdc1aInternalVariables::VFE);
        const auto EFD  = static_cast<size_t>(Esdc1aInternalVariables::EFD);

        const auto OMEGA = static_cast<size_t>(Esdc1aExternalVariables::OMEGA);
        const auto VREF  = static_cast<size_t>(Esdc1aExternalVariables::VREF);
        const auto VS    = static_cast<size_t>(Esdc1aExternalVariables::VS);
        const auto VUEL  = static_cast<size_t>(Esdc1aExternalVariables::VUEL);

        const ScalarT efdp = y[EFDP];
        const ScalarT vc   = y[VC];
        const ScalarT vr   = y[VR];
        const ScalarT vf   = y[VF];
        const ScalarT xll  = y[XLL];
        const ScalarT ev   = y[EV];
        const ScalarT vll  = y[VLL];
        const ScalarT vhv  = y[VHV];
        const ScalarT se   = y[SE];
        const ScalarT vfe  = y[VFE];
        const ScalarT efd  = y[EFD];

        const ScalarT efdp_dot = yp[EFDP];
        const ScalarT vc_dot   = yp[VC];
        const ScalarT vr_dot   = yp[VR];
        const ScalarT vf_dot   = yp[VF];
        const ScalarT xll_dot  = yp[XLL];

        const ScalarT omega = ws[OMEGA];
        const ScalarT vref  = ws[VREF];
        const ScalarT vs    = ws[VS];
        const ScalarT vuel  = ws[VUEL];

        const ScalarT ec                = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]);
        const ScalarT ev_target         = vref + vs + uel_on_ * vuel - vc - vf;
        const ScalarT vfe_target        = (Ke_ + se) * efdp;
        const ScalarT efdp_rate         = (vr - vfe) / Te_;
        const ScalarT limited_efdp_rate = awmin(efdp, efdp_rate, ZERO<RealT>);

        f[EFDP] = -efdp_dot + (ONE<RealT> - lim_on_) * efdp_rate
                  + lim_on_ * limited_efdp_rate;
        f[VC]  = -vc_dot + (ec - vc) / Tr_;
        f[VR]  = -vr_dot + Math::antiwindup(vr, -vr + Ka_ * vhv, Vrmin_, Vrmax_) / Ta_;
        f[VF]  = -vf_dot + (-vf + Kf_ * (vr - vfe) / Te_) / Tf1_;
        f[XLL] = -xll_dot + (ev - xll) / Tb_;
        f[EV]  = -ev + ev_target;
        f[VLL] = -vll + xll + (Tc_ / Tb_) * (ev - xll);
        f[VHV] = -vhv + uel_on_ * vll
                 + (ONE<RealT> - uel_on_) * Math::max(vll, vuel);
        f[SE]  = -se + SB_ * Math::qramp(efdp - SA_);
        f[VFE] = -vfe + vfe_target;
        f[EFD] = -efd + (ONE<RealT> + spd_on_ * omega) * efdp;

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Smooth anti-windup derivative above a fixed lower bound
       *
       * Passes the unconstrained rate above the bound, admits restoring
       * motion from below it, and smoothly blocks outward motion.
       *
       * @param[in] x State limited from below.
       * @param[in] f Unconstrained derivative of @p x.
       * @param[in] xmin Fixed lower bound on @p x.
       * @return Anti-windup-limited derivative.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline scalar_type
      Esdc1a<scalar_type, index_type>::awmin(
          const ScalarT x,
          const ScalarT f,
          const RealT   xmin)
      {
        const ScalarT above = Math::above(x, xmin);

        return (above + (ONE<RealT> - above) * Math::sigmoid(f)) * f;
      }

      /**
       * @brief Read the parameters out of the model data
       *
       * No parameter is required; every parameter keeps the default
       * documented in the model README when omitted. A non-numeric value, a
       * switch outside \f$\{0,1\}\f$, or a non-integer selector is counted and
       * reported by verify() rather than throwing. Integer JSON values are
       * accepted for real parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void Esdc1a<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        auto load_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            target = *real_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            target = static_cast<RealT>(*index_value);
          }
          else
          {
            Log::error() << "Esdc1a: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto load_switch = [&](auto key, bool& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* bool_value = std::get_if<bool>(&value))
          {
            target = *bool_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value);
                   index_value && (*index_value == 0 || *index_value == 1))
          {
            target = (*index_value == 1);
          }
          else if (const auto* real_value = std::get_if<RealT>(&value);
                   real_value && (*real_value == ZERO<RealT> || *real_value == ONE<RealT>) )
          {
            target = (*real_value == ONE<RealT>);
          }
          else
          {
            Log::error() << "Esdc1a: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        auto load_selector = [&](auto key, IdxT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            target = *index_value;
          }
          else if (const auto* real_value = std::get_if<RealT>(&value);
                   real_value && *real_value >= ZERO<RealT>
                   && *real_value == std::round(*real_value))
          {
            target = static_cast<IdxT>(std::round(*real_value));
          }
          else
          {
            Log::error() << "Esdc1a: parameter '" << name << "' must be an integer selector\n";
            ++parameter_error_count_;
          }
        };

        load_real(Params::Tr, Tr_, "Tr");
        load_real(Params::Ka, Ka_, "Ka");
        load_real(Params::Ta, Ta_, "Ta");
        load_real(Params::Tb, Tb_, "Tb");
        load_real(Params::Tc, Tc_, "Tc");
        load_real(Params::Vrmax, Vrmax_, "Vrmax");
        load_real(Params::Vrmin, Vrmin_, "Vrmin");
        load_real(Params::Ke, Ke_, "Ke");
        load_real(Params::Te, Te_, "Te");
        load_real(Params::Kf, Kf_, "Kf");
        load_real(Params::Tf1, Tf1_, "Tf1");
        load_switch(Params::Spdmlt, Spdmlt_, "Spdmlt");
        load_real(Params::E1, E1_, "E1");
        load_real(Params::Se1, Se1_, "Se1");
        load_real(Params::E2, E2_, "E2");
        load_real(Params::Se2, Se2_, "Se2");
        load_selector(Params::UEL, UEL_, "UEL");
        load_switch(Params::exclim, exclim_, "exclim");
        setDerivedParameters();
      }

      /**
       * @brief Bind the monitorable variables to their internal states
       *
       * Binds configured monitor keys to their corresponding internal
       * variables.
       */
      template <typename scalar_type, typename index_type>
      void Esdc1a<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::efd, [this]
                      { return y_.getData()[static_cast<size_t>(Esdc1aInternalVariables::EFD)]; });
        monitor_->set(Variable::vc, [this]
                      { return y_.getData()[static_cast<size_t>(Esdc1aInternalVariables::VC)]; });
        monitor_->set(Variable::vr, [this]
                      { return y_.getData()[static_cast<size_t>(Esdc1aInternalVariables::VR)]; });
        monitor_->set(Variable::vf, [this]
                      { return y_.getData()[static_cast<size_t>(Esdc1aInternalVariables::VF)]; });
        monitor_->set(Variable::se, [this]
                      { return y_.getData()[static_cast<size_t>(Esdc1aInternalVariables::SE)]; });
        monitor_->set(Variable::vfe, [this]
                      { return y_.getData()[static_cast<size_t>(Esdc1aInternalVariables::VFE)]; });
      }

      /**
       * @brief Resolve the parameter-derived constants and selector masks
       *
       * Raises the transducer, regulator, lead-lag, exciter, and feedback
       * lags to the well-posedness floor, fits the quadratic saturation
       * curve, and turns the three selectors into multiplicative masks. The
       * masks let the residual select signal routing without
       * parameter-dependent control flow, which keeps its structure fixed for
       * sparse automatic differentiation.
       */
      template <typename scalar_type, typename index_type>
      void Esdc1a<scalar_type, index_type>::setDerivedParameters()
      {
        // The lags are raised to the floor in place, so a negative value is
        // rejected here while the value as read is still available. verify()
        // reports the count.
        auto check_non_negative = [&](RealT value, const char* name)
        {
          if (value < ZERO<RealT>)
          {
            Log::error() << "Esdc1a: " << name << " must be non-negative\n";
            ++parameter_error_count_;
          }
        };

        check_non_negative(Tr_, "Tr");
        check_non_negative(Ta_, "Ta");
        check_non_negative(Tb_, "Tb");
        check_non_negative(Te_, "Te");
        check_non_negative(Tf1_, "Tf1");

        if (Tr_ < TIME_CONSTANT_MINIMUM || Ta_ < TIME_CONSTANT_MINIMUM
            || Tb_ < TIME_CONSTANT_MINIMUM || Te_ < TIME_CONSTANT_MINIMUM
            || Tf1_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "Esdc1a: Tr, Ta, Tb, Te, and Tf1 below "
                         << TIME_CONSTANT_MINIMUM
                         << " s are raised to that floor to keep the exciter lags well posed\n";
        }

        Tr_  = std::max(Tr_, TIME_CONSTANT_MINIMUM);
        Ta_  = std::max(Ta_, TIME_CONSTANT_MINIMUM);
        Tb_  = std::max(Tb_, TIME_CONSTANT_MINIMUM);
        Te_  = std::max(Te_, TIME_CONSTANT_MINIMUM);
        Tf1_ = std::max(Tf1_, TIME_CONSTANT_MINIMUM);

        spd_on_ = ZERO<RealT>;
        if (Spdmlt_)
        {
          spd_on_ = ONE<RealT>;
        }

        uel_on_ = ZERO<RealT>;
        if (UEL_ >= static_cast<IdxT>(2))
        {
          uel_on_ = ONE<RealT>;
        }

        lim_on_ = ZERO<RealT>;
        if (exclim_)
        {
          lim_on_ = ONE<RealT>;
        }

        // A disabled or inconsistent saturation curve keeps the zero fit so
        // the coefficients stay finite; verify() reports inconsistent data.
        const bool saturation_enabled = !(Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>);
        const bool saturation_points_are_ordered =
            (E2_ > E1_ && Se2_ > Se1_)
            || (E2_ < E1_ && Se2_ < Se1_);
        const bool saturation_consistent =
            E1_ > ZERO<RealT> && E2_ > ZERO<RealT>
            && Se1_ > ZERO<RealT> && Se2_ > ZERO<RealT>
            && saturation_points_are_ordered;
        if (!saturation_enabled || !saturation_consistent)
        {
          SA_ = ZERO<RealT>;
          SB_ = ZERO<RealT>;
          return;
        }

        const RealT C = std::sqrt(Se2_ / Se1_);
        SA_           = (C * E1_ - E2_) / (C - ONE<RealT>);
        SB_           = Se1_ / ((E1_ - SA_) * (E1_ - SA_));
      }

      /**
       * @brief Invert the smooth CommonMath ramp
       *
       * Initialization seeds the inactive high-value gate with the gate
       * *input*, so the residual reproduces the requested output through the
       * same smooth ramp it evaluates. Beyond the softplus width the smooth
       * ramp is the identity to double precision, so the output is returned
       * unchanged there.
       *
       * @param[in] ramp_output Strictly positive requested ramp output.
       * @return The input the smooth ramp maps to the requested output.
       *
       * @pre @p ramp_output is finite and strictly positive.
       * @warning This function contains conditional branching and may be used
       *          during initialization, but not during residual or Jacobian
       *          evaluation.
       */
      template <typename scalar_type, typename index_type>
      typename Esdc1a<scalar_type, index_type>::RealT
      Esdc1a<scalar_type, index_type>::inverseRamp(RealT ramp_output) const
      {
        static constexpr RealT SOFTPLUS_WIDTH = static_cast<RealT>(50.0);

        const RealT scaled_output = Math::MU<RealT> * ramp_output;
        if (scaled_output > SOFTPLUS_WIDTH)
        {
          return ramp_output;
        }
        return std::log(std::expm1(scaled_output)) / Math::MU<RealT>;
      }

      /**
       * @brief Access the terminal-bus real voltage component.
       *
       * @return Reference to the bus variable.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Esdc1a<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      /**
       * @brief Access the terminal-bus imaginary voltage component.
       *
       * @return Reference to the bus variable.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Esdc1a<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

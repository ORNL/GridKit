/**
 * @file HygovImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the HYGOV governor model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct a HYGOV governor without parameters
       *
       * The model is sized but left unconfigured. Every parameter keeps its
       * documented default, the required power base is absent, and no monitor
       * is created, so verify() reports configuration errors until the data
       * constructor is used instead.
       */
      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::Hygov()
      {
        size_ = static_cast<IdxT>(HygovIdx::MAXIMUM);
      }

      /**
       * @brief Construct a HYGOV governor from model data
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::Hygov(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(HygovIdx::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::~Hygov()
      {
      }

      /**
       * @brief Resolve the parameter-derived constants
       *
       * Raises each governor lag to the well-posedness floor, sizes the
       * component power base, and derives the speed lead-lag gain from the
       * floored denominator so the residual keeps a fixed structure for
       * sparse automatic differentiation.
       */
      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::setDerivedParameters()
      {
        // The lags are raised to the floor in place, so a negative value is
        // rejected here while the value as read is still available. verify()
        // reports the count.
        auto check_non_negative = [&](RealT value, const char* name)
        {
          if (value < ZERO<RealT>)
          {
            Log::error() << "Hygov: " << name << " must be non-negative\n";
            ++parameter_error_count_;
          }
        };

        check_non_negative(Tr_, "Tr");
        check_non_negative(Tf_, "Tf");
        check_non_negative(Tg_, "Tg");
        check_non_negative(Tw_, "Tw");
        check_non_negative(Tnp_, "Tnp");

        if (Tr_ < TIME_CONSTANT_MINIMUM || Tf_ < TIME_CONSTANT_MINIMUM
            || Tg_ < TIME_CONSTANT_MINIMUM || Tw_ < TIME_CONSTANT_MINIMUM
            || Tnp_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "Hygov: Tr, Tf, Tg, Tw, and Tnp below "
                         << TIME_CONSTANT_MINIMUM
                         << " s are raised to that floor to keep the governor lags well posed\n";
        }

        Tr_  = std::max(Tr_, TIME_CONSTANT_MINIMUM);
        Tf_  = std::max(Tf_, TIME_CONSTANT_MINIMUM);
        Tg_  = std::max(Tg_, TIME_CONSTANT_MINIMUM);
        Tw_  = std::max(Tw_, TIME_CONSTANT_MINIMUM);
        Tnp_ = std::max(Tnp_, TIME_CONSTANT_MINIMUM);

        va_component_base_ = Trate_ * static_cast<RealT>(1.0e6);

        leadlag_gain_ = Tn_ / Tnp_;
      }

      /**
       * @brief Evaluate the nonlinear gate-to-power curve
       *
       * Sums the five smooth CommonMath linear segments spanned by the
       * `Gv`/`Pgv` points, so the same fixed expression serves the residual
       * and both scalar instantiations.
       *
       * @param[in] gate Gate position.
       * @return Turbine power at nominal head.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline scalar_type
      Hygov<scalar_type, index_type>::gatePower(scalar_type gate) const
      {
        return ScalarT{Pgv_[0]}
               + Math::linseg(gate, Gv_[0], Gv_[1], Pgv_[1] - Pgv_[0])
               + Math::linseg(gate, Gv_[1], Gv_[2], Pgv_[2] - Pgv_[1])
               + Math::linseg(gate, Gv_[2], Gv_[3], Pgv_[3] - Pgv_[2])
               + Math::linseg(gate, Gv_[3], Gv_[4], Pgv_[4] - Pgv_[3])
               + Math::linseg(gate, Gv_[4], Gv_[5], Pgv_[5] - Pgv_[4]);
      }

      /**
       * @brief Solve the steady gate position for a seeded mechanical power
       *
       * At the steady state the head rests at the dam head, the flow rides
       * the gate curve, and the turbine power less the speed-damping loss
       * reproduces the seed:
       * @f[
       *   A_t H_0 \left(\sqrt{H_0}\,N_{\mathrm{GV}}(g) - q_{\mathrm{NL}}\right)
       *   - D_{\mathrm{turb}}\,\omega_0\, g = P_{\mathrm{m},0}.
       * @f]
       * The curve is linear on each rising segment, so the equation is solved
       * segment by segment and the lowest admissible gate is selected. Flat
       * segments carry no power information and are skipped.
       *
       * @param[in] pmech Seeded mechanical power on the component base.
       * @param[in] omega Initial machine speed deviation.
       * @return The gate position, or a quiet NaN when no rising segment
       *         reproduces the seed.
       */
      template <typename scalar_type, typename index_type>
      typename Hygov<scalar_type, index_type>::RealT
      Hygov<scalar_type, index_type>::solveInitialGate(RealT pmech, RealT omega) const
      {
        const RealT h0      = Hdam_;
        const RealT gain    = At_ * h0 * std::sqrt(h0);
        const RealT damping = Dturb_ * omega;
        const RealT target  = pmech + At_ * h0 * Qnl_;

        if (std::abs(gain * Pgv_[0] - damping * Gv_[0] - target) <= INITIALIZATION_TOLERANCE)
        {
          return Gv_[0];
        }

        for (size_t i = 0; i < 5; ++i)
        {
          if (Pgv_[i + 1] <= Pgv_[i])
          {
            continue;
          }

          const RealT slope       = (Pgv_[i + 1] - Pgv_[i]) / (Gv_[i + 1] - Gv_[i]);
          const RealT denominator = gain * slope - damping;
          if (std::abs(denominator) <= INITIALIZATION_TOLERANCE)
          {
            continue;
          }

          const RealT gate = (target - gain * (Pgv_[i] - slope * Gv_[i])) / denominator;
          if (Gv_[i] - INITIALIZATION_TOLERANCE <= gate
              && gate <= Gv_[i + 1] + INITIALIZATION_TOLERANCE)
          {
            return gate;
          }
        }

        return std::numeric_limits<RealT>::quiet_NaN();
      }

      /**
       * @brief Convert a system-base power to HYGOV component base
       *
       * @param[in] value Quantity on the system base.
       * @return The same quantity on the component base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Hygov<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * va_system_base_ / va_component_base_;
      }

      /**
       * @brief Convert a component-base power to the system base
       *
       * @param[in] value Quantity on the component base.
       * @return The same quantity on the system base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Hygov<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      /**
       * @brief Read the parameters out of the model data
       *
       * Only the turbine-rating power base is required; every other parameter
       * keeps the default documented in the model README when omitted. A
       * missing required key or a non-numeric value is counted and reported
       * by verify() rather than throwing. Integer JSON values are accepted
       * for real parameters. All-zero `Gv` and `Pgv` source points select the
       * identity gate curve.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
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
            Log::error() << "Hygov: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        if (!data.parameters.contains(Params::Trate))
        {
          Log::error() << "Hygov: missing required parameter 'Trate'\n";
          ++parameter_error_count_;
        }
        load_real(Params::Trate, Trate_, "Trate");
        load_real(Params::Rperm, Rperm_, "Rperm");
        load_real(Params::Rtemp, Rtemp_, "Rtemp");
        load_real(Params::Tr, Tr_, "Tr");
        load_real(Params::Tf, Tf_, "Tf");
        load_real(Params::Tg, Tg_, "Tg");
        load_real(Params::Velm, Velm_, "Velm");
        load_real(Params::Gmax, Gmax_, "Gmax");
        load_real(Params::Gmin, Gmin_, "Gmin");
        load_real(Params::Tw, Tw_, "Tw");
        load_real(Params::At, At_, "At");
        load_real(Params::Dturb, Dturb_, "Dturb");
        load_real(Params::Qnl, Qnl_, "Qnl");
        load_real(Params::Tn, Tn_, "Tn");
        load_real(Params::Tnp, Tnp_, "Tnp");
        load_real(Params::db1, db1_, "db1");
        load_real(Params::db2, db2_, "db2");
        load_real(Params::Hdam, Hdam_, "Hdam");
        load_real(Params::Gv0, Gv_[0], "Gv0");
        load_real(Params::Gv1, Gv_[1], "Gv1");
        load_real(Params::Gv2, Gv_[2], "Gv2");
        load_real(Params::Gv3, Gv_[3], "Gv3");
        load_real(Params::Gv4, Gv_[4], "Gv4");
        load_real(Params::Gv5, Gv_[5], "Gv5");
        load_real(Params::Pgv0, Pgv_[0], "Pgv0");
        load_real(Params::Pgv1, Pgv_[1], "Pgv1");
        load_real(Params::Pgv2, Pgv_[2], "Pgv2");
        load_real(Params::Pgv3, Pgv_[3], "Pgv3");
        load_real(Params::Pgv4, Pgv_[4], "Pgv4");
        load_real(Params::Pgv5, Pgv_[5], "Pgv5");

        const bool source_default_curve =
            std::all_of(Gv_.begin(), Gv_.end(), [](RealT value)
                        { return value == ZERO<RealT>; })
            && std::all_of(Pgv_.begin(), Pgv_.end(), [](RealT value)
                           { return value == ZERO<RealT>; });
        if (source_default_curve)
        {
          Gv_  = {ZERO<RealT>,
                  static_cast<RealT>(0.2),
                  static_cast<RealT>(0.4),
                  static_cast<RealT>(0.6),
                  static_cast<RealT>(0.8),
                  ONE<RealT>};
          Pgv_ = Gv_;
        }

        setDerivedParameters();
      }

      /**
       * @brief Access the monitor
       *
       * @return Monitor for this model, or nullptr when the model was
       *         constructed without data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Hygov<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Bind the monitorable variables to their internal states
       *
       * The mechanical-power output is published on the system base and the
       * remaining outputs on the component base, as documented in the model
       * README.
       */
      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::initializeMonitor()
      {
        using I        = HygovIdx;
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::pmech, [this]
                      { return y_.getData()[I::PMECH]; });
        monitor_->set(Variable::filter, [this]
                      { return y_.getData()[I::XF]; });
        monitor_->set(Variable::desiredgate, [this]
                      { return y_.getData()[I::C]; });
        monitor_->set(Variable::gate, [this]
                      { return y_.getData()[I::G]; });
        monitor_->set(Variable::flow, [this]
                      { return y_.getData()[I::Q]; });
        monitor_->set(Variable::head, [this]
                      { return y_.getData()[I::H]; });
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the model vectors and wire the mechanical-power output
       *
       * Sizes the state, residual, and signal-interface buffers, seeds the
       * identity index maps, and points the assigned `pmech` node at the
       * internal state it publishes. That node aliases HYGOV storage from
       * here on, which is how initialize() reads the seed the machine wrote.
       * HYGOV attaches to no bus, so the bus-interface buffer stays empty.
       * Repeated calls reuse the allocated vectors.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::allocate()
      {
        using I = HygovIdx;
        using E = HygovExt;

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.clear();

        auto signal_size = E::MAXIMUM;
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<HygovInternalVariables::PMECH>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<HygovInternalVariables::PMECH>()->set(
              &y[I::PMECH],
              &(this->getVariableIndex(static_cast<IdxT>(I::PMECH))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the HYGOV configuration
       *
       * Checks parameter-loading errors, static parameter relationships, the
       * gate-curve monotonicity, and attached external signals. Seed
       * feasibility is operating-point dependent and is checked by
       * initialize().
       *
       * @return int Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Hygov: " << message << '\n';
            ret += 1;
          }
        };

        check(Trate_ > ZERO<RealT>, "Trate must be positive");
        check(Rtemp_ != ZERO<RealT>, "Rtemp must be nonzero");
        check(Tn_ >= ZERO<RealT>, "Tn must be non-negative");
        check(Velm_ >= ZERO<RealT>, "Velm must be non-negative");
        check(Gmin_ <= Gmax_, "Gmin must be less than or equal to Gmax");
        check(At_ > ZERO<RealT>, "At must be positive");
        check(Dturb_ >= ZERO<RealT>, "Dturb must be non-negative");
        check(db1_ >= ZERO<RealT>, "db1 must be non-negative");
        check(Hdam_ > ZERO<RealT>, "Hdam must be positive");

        for (size_t i = 1; i < Gv_.size(); ++i)
        {
          check(Gv_[i - 1] < Gv_[i], "Gv points must be strictly increasing");
          check(Pgv_[i - 1] <= Pgv_[i], "Pgv points must be non-decreasing");
        }

        // An attached port must resolve to readable signal storage. The
        // enumerator is a template argument, so each port names itself once.
        auto check_attached_signal =
            [&]<HygovExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "Hygov: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_attached_signal.template operator()<HygovExternalVariables::OMEGA>("speed");
        check_attached_signal.template operator()<HygovExternalVariables::PREF>("pref");
        check_attached_signal.template operator()<HygovExternalVariables::PAUX>("paux");

        return ret;
      }

      /**
       * @brief Initialize HYGOV from the seeded mechanical-power port
       *
       * Reads the assigned system-base `pmech` node and the attached speed
       * and auxiliary-power inputs, solves the component-base steady state
       * that preserves the seed, and publishes the resolved load reference to
       * an attached `pref` signal. All operating-point checks are completed
       * before model or signal storage is modified.
       *
       * @pre allocate() has completed.
       * @pre The machine model has seeded the assigned `pmech` node.
       *
       * @return int 0 on success; nonzero when the configuration is invalid,
       *             no rising segment of the gate curve reproduces the seeded
       *             power, or the resulting gate is outside Gmin/Gmax.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::initialize()
      {
        using I = HygovIdx;

        if (verify() > 0)
        {
          Log::error() << "Hygov: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y = y_.getData();

        // The assigned pmech node aliases this entry after allocate(). Its
        // system-base seed remains untouched throughout initialization.
        const ScalarT pmech0 = toComponentBase(y[I::PMECH]);

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<HygovExternalVariables::OMEGA>())
        {
          omega0 = signals_.template readExternalVariable<HygovExternalVariables::OMEGA>();
        }

        ScalarT paux0_system{ZERO<RealT>};
        if (signals_.template isAttached<HygovExternalVariables::PAUX>())
        {
          paux0_system = signals_.template readExternalVariable<HygovExternalVariables::PAUX>();
        }
        const ScalarT paux0 = toComponentBase(paux0_system);

        const RealT gate0 = solveInitialGate(static_cast<RealT>(pmech0),
                                             static_cast<RealT>(omega0));
        if (std::isnan(gate0))
        {
          Log::error() << "Hygov: initial mechanical power is outside the invertible gate curve\n";
          return 1;
        }
        if (gate0 < Gmin_ - INITIALIZATION_TOLERANCE
            || gate0 > Gmax_ + INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Hygov: initialized gate is outside Gmin/Gmax\n";
          return 1;
        }

        const ScalarT h0       = static_cast<ScalarT>(Hdam_);
        const ScalarT pgv0     = gatePower(static_cast<ScalarT>(gate0));
        const ScalarT q0       = std::sqrt(Hdam_) * pgv0;
        const ScalarT omegadb0 = Math::deadband1(omega0, -db1_, db1_);
        const ScalarT xn0      = omegadb0;
        const ScalarT yomega0  = xn0 + leadlag_gain_ * (omegadb0 - xn0);
        const ScalarT pref0    = toSystemBase(yomega0 + Rperm_ * gate0 - paux0);

        y[I::XN]      = xn0;
        y[I::XF]      = ZERO<RealT>;
        y[I::C]       = gate0;
        y[I::G]       = gate0;
        y[I::Q]       = q0;
        y[I::OMEGADB] = omegadb0;
        y[I::EF]      = ZERO<RealT>;
        y[I::FC]      = ZERO<RealT>;
        y[I::RC]      = ZERO<RealT>;
        y[I::PGV]     = pgv0;
        y[I::H]       = h0;

        pref_set_ = pref0;
        paux_set_ = paux0_system;

        if (signals_.template isAttached<HygovExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<HygovExternalVariables::PREF>(pref_set_);
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief Identify the differential variables
       *
       * The speed lead-lag state, the governor error filter, the desired
       * gate, the gate servo, and the turbine flow carry derivatives; every
       * other internal variable is algebraic.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::tagDifferentiable()
      {
        using I = HygovIdx;

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[I::XN] = true;
        tag_[I::XF] = true;
        tag_[I::C]  = true;
        tag_[I::G]  = true;
        tag_[I::Q]  = true;
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * All HYGOV variables are per-unit speeds, gates, flows, heads, and
       * powers of the same order, so they share the relative tolerance as
       * their absolute floor.
       *
       * @param[in] rel_tol Solver relative tolerance.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Internal residual
       *
       * Evaluates the five governor states and the seven algebraic rows
       * documented in the model README. The body is kept free of branches
       * and loops so that sparse automatic differentiation resolves a fixed
       * structure; the gate curve enters as a fixed sum of smooth linear
       * segments.
       *
       * @param[in] y Internal variables.
       * @param[in] yp Internal variable derivatives.
       * @param[in] wb Bus voltage components; unused, HYGOV attaches to no bus.
       * @param[in] ws External signal values on system base.
       * @param[out] f Internal residuals.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Hygov<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        using I = HygovIdx;
        using E = HygovExt;

        const ScalarT xn      = y[I::XN];
        const ScalarT xf      = y[I::XF];
        const ScalarT c       = y[I::C];
        const ScalarT g       = y[I::G];
        const ScalarT q       = y[I::Q];
        const ScalarT omegadb = y[I::OMEGADB];
        const ScalarT ef      = y[I::EF];
        const ScalarT fc      = y[I::FC];
        const ScalarT rc      = y[I::RC];
        const ScalarT pgv     = y[I::PGV];
        const ScalarT head    = y[I::H];
        const ScalarT pmech   = y[I::PMECH];

        const ScalarT xn_dot = yp[I::XN];
        const ScalarT xf_dot = yp[I::XF];
        const ScalarT c_dot  = yp[I::C];
        const ScalarT g_dot  = yp[I::G];
        const ScalarT q_dot  = yp[I::Q];

        const ScalarT omega = ws[E::OMEGA];
        const ScalarT pref  = toComponentBase(ws[E::PREF]);
        const ScalarT paux  = toComponentBase(ws[E::PAUX]);

        const ScalarT yomega = xn + leadlag_gain_ * (omegadb - xn);

        f[I::XN]      = -xn_dot + (omegadb - xn) / Tnp_;
        f[I::XF]      = -xf_dot + (ef - xf) / Tf_;
        f[I::C]       = -c_dot + Math::antiwindup(c, rc, Gmin_, Gmax_);
        f[I::G]       = -g_dot + (c - g) / Tg_;
        f[I::Q]       = -q_dot + (Hdam_ - head) / Tw_;
        f[I::OMEGADB] = -omegadb + Math::deadband1(omega, -db1_, db1_);
        f[I::EF]      = -ef + pref + paux - yomega - Rperm_ * c;
        f[I::FC]      = -fc + (xf / Tr_ + (ef - xf) / Tf_) / Rtemp_;
        f[I::RC]      = -rc + Math::clamp(fc, -Velm_, Velm_);
        f[I::PGV]     = -pgv + gatePower(g);
        f[I::H]       = -q * q + head * pgv * pgv;
        f[I::PMECH]   = -toComponentBase(pmech) + At_ * head * (q - Qnl_) - Dturb_ * omega * g;

        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       * Refreshes the signal interface buffers and evaluates the internal
       * residual. HYGOV attaches to no bus, so there is no bus interface to
       * refresh. An unattached reference or auxiliary port falls back to the
       * value latched by initialize(); an unattached speed port reads zero
       * deviation.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::evaluateResidual()
      {
        using E = HygovExt;

        ws_[E::OMEGA] = ZERO<RealT>;
        ws_[E::PREF]  = pref_set_;
        ws_[E::PAUX]  = paux_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<HygovExternalVariables::OMEGA>())
        {
          ws_[E::OMEGA] = signals_.template readExternalVariable<HygovExternalVariables::OMEGA>();
          ws_indices_[E::OMEGA] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<HygovExternalVariables::PREF>())
        {
          ws_[E::PREF] = signals_.template readExternalVariable<HygovExternalVariables::PREF>();
          ws_indices_[E::PREF] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::PREF>();
        }
        if (signals_.template isAttached<HygovExternalVariables::PAUX>())
        {
          ws_[E::PAUX] = signals_.template readExternalVariable<HygovExternalVariables::PAUX>();
          ws_indices_[E::PAUX] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::PAUX>();
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

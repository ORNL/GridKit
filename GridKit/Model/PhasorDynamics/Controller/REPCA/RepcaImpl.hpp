/**
 * @file RepcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REPCA phasor-dynamics plant-control model.
 */

#pragma once

#include <algorithm>
#include <mutex>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Controller
    {
      /// Logger used for REPCA diagnostics.
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct REPCA with its documented parameter defaults
       *
       * The regulated bus is retained, the model is sized, and every
       * parameter keeps its documented default. No monitor or signal
       * connection is created.
       *
       * @param[in] bus Regulated bus measured by the controller.
       */
      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(RepcaInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      /**
       * @brief Construct REPCA from model data
       *
       * @param[in] bus Regulated bus measured by the controller.
       * @param[in] data Model parameters and monitor selections.
       */
      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(RepcaInternalVariables::MAXIMUM);
      }

      /**
       * @brief Destroy the plant controller and its optional variable monitor.
       */
      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::~Repca()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate model vectors and wire assigned command outputs
       *
       * Sizes the state, residual, bus, and signal-interface buffers, initializes
       * identity index maps, and points assigned `qext` and `pext` nodes at
       * the internal system-base states that REPCA publishes. Repeated
       * allocation reuses the existing model vectors and signal links.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::allocate()
      {
        const auto QEXT = static_cast<size_t>(RepcaInternalVariables::QEXT);
        const auto PEXT = static_cast<size_t>(RepcaInternalVariables::PEXT);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        const auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(2, ScalarT{0});

        const auto signal_size = static_cast<size_t>(RepcaExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        auto* y = y_.getData();

        if (signals_.template isAssigned<RepcaInternalVariables::QEXT>())
        {
          signals_.template getSignalNode<RepcaInternalVariables::QEXT>()->set(
              &y[QEXT],
              &(this->getVariableIndex(static_cast<IdxT>(QEXT))));
        }

        if (signals_.template isAssigned<RepcaInternalVariables::PEXT>())
        {
          signals_.template getSignalNode<RepcaInternalVariables::PEXT>()->set(
              &y[PEXT],
              &(this->getVariableIndex(static_cast<IdxT>(PEXT))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the REPCA configuration
       *
       * Checks parameter-loading errors, static parameter relationships,
       * system/component bases and both conversion ratios, the regulated
       * bus, required measurement signals, and attached optional reference
       * signals. Command-output assignment is optional.
       *
       * @return Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Repca: " << message << '\n';
            ret += 1;
          }
        };

        check(bus_ != nullptr, "regulated bus is required");

        const bool valid_component_base = std::isfinite(mva_base_)
                                          && mva_base_ > ZERO<RealT>
                                          && std::isfinite(va_component_base_)
                                          && va_component_base_ > ZERO<RealT>;
        const bool valid_system_base = std::isfinite(va_system_base_)
                                       && va_system_base_ > ZERO<RealT>;
        check(valid_component_base,
              "mva must define a finite positive component power base");
        check(valid_system_base, "system power base must be finite and positive");
        if (valid_component_base && valid_system_base)
        {
          const RealT system_to_component = va_system_base_ / va_component_base_;
          const RealT component_to_system = va_component_base_ / va_system_base_;
          check(std::isfinite(system_to_component)
                    && system_to_component > ZERO<RealT>
                    && std::isfinite(component_to_system)
                    && component_to_system > ZERO<RealT>,
                "system/component power-base conversion ratios must be finite and positive");
        }

        check(dbdlow_ <= ZERO<RealT> && ZERO<RealT> <= dbdupper_,
              "dbdlow <= 0 <= dbdupper is required");
        check(emin_ <= ZERO<RealT> && ZERO<RealT> <= emax_,
              "emin <= 0 <= emax is required");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(fdbd1_ <= ZERO<RealT> && ZERO<RealT> <= fdbd2_,
              "fdbd1 <= 0 <= fdbd2 is required");
        check(Ddn_ >= ZERO<RealT>, "Ddn must be non-negative");
        check(Dup_ >= ZERO<RealT>, "Dup must be non-negative");
        check(femin_ <= ZERO<RealT> && ZERO<RealT> <= femax_,
              "femin <= 0 <= femax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");

        auto check_required_signal = [&]<RepcaExternalVariables variable>(const char* name)
        {
          if (!signals_.template isAttached<variable>())
          {
            Log::error() << "Repca: " << name << " signal is required\n";
            ret += 1;
          }
          else if (!signals_.template isLinked<variable>())
          {
            Log::error() << "Repca: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_required_signal.template operator()<RepcaExternalVariables::IR>("ir");
        check_required_signal.template operator()<RepcaExternalVariables::II>("ii");
        check_required_signal.template operator()<RepcaExternalVariables::P>("p");
        check_required_signal.template operator()<RepcaExternalVariables::Q>("q");

        auto check_optional_signal = [&]<RepcaExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "Repca: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_optional_signal.template operator()<RepcaExternalVariables::VREF>("vref");
        check_optional_signal.template operator()<RepcaExternalVariables::PREF>("pref");
        check_optional_signal.template operator()<RepcaExternalVariables::QREF>("qref");
        check_optional_signal.template operator()<RepcaExternalVariables::FREQ>("freq");
        check_optional_signal.template operator()<RepcaExternalVariables::FREQREF>("freqref");

        return ret;
      }

      /**
       * @brief Initialize REPCA from the initial plant commands
       *
       * Reads the required bus and measurement ports, preserves the initial
       * system-base `qext` value and, when frequency control is enabled, `pext`,
       * reconstructs the steady controller state, and publishes the resolved
       * optional references.
       *
       * @pre allocate() has completed.
       * @pre verify() reports a valid parameter and port configuration.
       * @pre `qext` and, when frequency control is enabled, `pext` contain the
       *      initial plant commands.
       *
       * @post On failure, states, derivatives, limits, latches, and signal
       *       values are unchanged.
       *
       * @return 0 on success; nonzero when allocation, configuration, initial-
       *         value, candidate, or steady-state checks fail.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::initialize()
      {
        const auto VMEAS  = static_cast<size_t>(RepcaInternalVariables::VMEAS);
        const auto QMEAS  = static_cast<size_t>(RepcaInternalVariables::QMEAS);
        const auto XQPI   = static_cast<size_t>(RepcaInternalVariables::XQPI);
        const auto XQLAG  = static_cast<size_t>(RepcaInternalVariables::XQLAG);
        const auto PMEAS  = static_cast<size_t>(RepcaInternalVariables::PMEAS);
        const auto XPPI   = static_cast<size_t>(RepcaInternalVariables::XPPI);
        const auto PREF   = static_cast<size_t>(RepcaInternalVariables::PREF);
        const auto V      = static_cast<size_t>(RepcaInternalVariables::V);
        const auto VLDC   = static_cast<size_t>(RepcaInternalVariables::VLDC);
        const auto VDROOP = static_cast<size_t>(RepcaInternalVariables::VDROOP);
        const auto VCTRL  = static_cast<size_t>(RepcaInternalVariables::VCTRL);
        const auto SFRZ   = static_cast<size_t>(RepcaInternalVariables::SFRZ);
        const auto ERQ    = static_cast<size_t>(RepcaInternalVariables::ERQ);
        const auto ERQDB  = static_cast<size_t>(RepcaInternalVariables::ERQDB);
        const auto ERQLIM = static_cast<size_t>(RepcaInternalVariables::ERQLIM);
        const auto QPI    = static_cast<size_t>(RepcaInternalVariables::QPI);
        const auto QEXT   = static_cast<size_t>(RepcaInternalVariables::QEXT);
        const auto EF     = static_cast<size_t>(RepcaInternalVariables::EF);
        const auto EP     = static_cast<size_t>(RepcaInternalVariables::EP);
        const auto EPLIM  = static_cast<size_t>(RepcaInternalVariables::EPLIM);
        const auto PPI    = static_cast<size_t>(RepcaInternalVariables::PPI);
        const auto PEXT   = static_cast<size_t>(RepcaInternalVariables::PEXT);

        if (!allocated_)
        {
          Log::error() << "Repca: allocate must complete before initialize\n";
          return 1;
        }

        if (verify() > 0)
        {
          Log::error() << "Repca: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y = y_.getData();

        const ScalarT qext0_system = y[QEXT];
        const ScalarT pext0_system = y[PEXT];
        const ScalarT qext0        = toComponentBase(qext0_system);
        const ScalarT pext0        = toComponentBase(pext0_system);

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();
        const ScalarT ir_system =
            signals_.template readExternalVariable<RepcaExternalVariables::IR>();
        const ScalarT ii_system =
            signals_.template readExternalVariable<RepcaExternalVariables::II>();
        const ScalarT p_system =
            signals_.template readExternalVariable<RepcaExternalVariables::P>();
        const ScalarT q_system =
            signals_.template readExternalVariable<RepcaExternalVariables::Q>();
        const ScalarT ir   = toComponentBase(ir_system);
        const ScalarT ii   = toComponentBase(ii_system);
        ScalarT       freq = static_cast<ScalarT>(ONE<RealT>);
        if (signals_.template isAttached<RepcaExternalVariables::FREQ>())
        {
          freq = signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();
        }

        auto is_finite = [](ScalarT value)
        {
          return std::isfinite(static_cast<RealT>(value));
        };
        if (!is_finite(vr) || !is_finite(vi) || !is_finite(ir_system)
            || !is_finite(ii_system) || !is_finite(p_system)
            || !is_finite(q_system) || !is_finite(freq) || !is_finite(qext0)
            || (Freqflag_ && !is_finite(pext0)))
        {
          Log::error() << "Repca: initial bus, signal, and command values must be finite\n";
          return 1;
        }

        const ScalarT p = toComponentBase(p_system);
        const ScalarT q = toComponentBase(q_system);

        const ScalarT vldc_r = vr - Rc_ * ir + Xc_ * ii;
        const ScalarT vldc_i = vi - Rc_ * ii - Xc_ * ir;

        const ScalarT v0      = std::sqrt(vr * vr + vi * vi);
        const ScalarT vldc0   = std::sqrt(vldc_r * vldc_r + vldc_i * vldc_i);
        const ScalarT vdroop0 = v0 + Kc_ * q;
        const ScalarT vctrl0  = vcomp_on_ * vldc0 + vcomp_off_ * vdroop0;
        const ScalarT vmeas0  = vctrl0;
        const ScalarT qmeas0  = q;
        const ScalarT pmeas0  = p;
        const ScalarT sfrz0   = Math::above(v0, Vfrz_);

        const ScalarT zero = static_cast<ScalarT>(ZERO<RealT>);
        ScalarT       erq0{};
        ScalarT       erqdb0{};
        if (!invertClamp(zero, emin_, emax_, erqdb0)
            || !invertDeadband(erqdb0, dbdlow_, dbdupper_, erq0))
        {
          Log::error() << "Repca: reactive error blocks have no finite steady input\n";
          return 1;
        }
        const ScalarT erqlim0 = zero;
        const ScalarT qpi0    = qext0;
        const ScalarT xqlag0  = qpi0;

        const RealT qmin = std::min(Qmin_, static_cast<RealT>(qpi0));
        const RealT qmax = std::max(Qmax_, static_cast<RealT>(qpi0));

        ScalarT qpi_input0{};
        if (!invertClamp(qpi0, qmin, qmax, qpi_input0))
        {
          Log::error() << "Repca: reactive-power PI limiter has no finite steady input\n";
          return 1;
        }
        const ScalarT xqpi0      = qpi_input0 - Kp_ * erqlim0;
        const ScalarT q_aw_rate0 = Math::antiwindup(qpi0, Ki_ * erqlim0, qmin, qmax);
        const ScalarT xqpi_rate0 = sfrz0 * q_aw_rate0;
        if (!is_finite(q_aw_rate0) || !is_finite(xqpi_rate0)
            || std::abs(static_cast<RealT>(xqpi_rate0)) > INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Repca: reactive-power PI state rate is nonzero at initialization\n";
          return 1;
        }

        ScalarT frequency_error0{};
        ScalarT ep0{};
        if (!invertDeadband(zero, fdbd1_, fdbd2_, frequency_error0)
            || !invertClamp(zero, femin_, femax_, ep0))
        {
          Log::error() << "Repca: active error blocks have no finite steady input\n";
          return 1;
        }
        const ScalarT ef0    = zero;
        const ScalarT pfreq0 = droop(ef0, Ddn_, Dup_);
        const ScalarT eplim0 = zero;
        const ScalarT pref0  = Freqflag_ ? pext0 : pmeas0;
        const ScalarT ppi0   = pref0;

        const RealT pmin = std::min(Pmin_, static_cast<RealT>(ppi0));
        const RealT pmax = std::max(Pmax_, static_cast<RealT>(ppi0));

        ScalarT ppi_input0{};
        if (!invertClamp(ppi0, pmin, pmax, ppi_input0))
        {
          Log::error() << "Repca: active-power PI limiter has no finite steady input\n";
          return 1;
        }
        const ScalarT xppi0      = ppi_input0 - Kpg_ * eplim0;
        const ScalarT p_aw_rate0 = Math::antiwindup(ppi0, Kig_ * eplim0, pmin, pmax);
        if (!is_finite(p_aw_rate0)
            || std::abs(static_cast<RealT>(p_aw_rate0)) > INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Repca: active-power PI antiwindup rate is nonzero at initialization\n";
          return 1;
        }

        const ScalarT pext_output0 = Freqflag_ ? pext0_system : zero;
        const ScalarT freqref0     = freq + frequency_error0;
        ScalarT       vref0        = vmeas0;
        ScalarT       qref0_system = q_system;
        if (RefFlag_)
        {
          vref0 += erq0;
        }
        else
        {
          qref0_system = toSystemBase(qmeas0 + erq0);
        }
        const ScalarT pref0_system = p_system + toSystemBase(ep0 - pfreq0);

        const bool candidates_are_finite =
            is_finite(qext0_system)
            && (!Freqflag_ || is_finite(pext0_system))
            && is_finite(qext0)
            && (!Freqflag_ || is_finite(pext0))
            && is_finite(p)
            && is_finite(q)
            && is_finite(vldc_r)
            && is_finite(vldc_i)
            && is_finite(v0)
            && is_finite(vldc0)
            && is_finite(vdroop0)
            && is_finite(vctrl0)
            && is_finite(vmeas0)
            && is_finite(qmeas0)
            && is_finite(pmeas0)
            && is_finite(sfrz0)
            && is_finite(erq0)
            && is_finite(erqdb0)
            && is_finite(erqlim0)
            && is_finite(qpi0)
            && is_finite(xqlag0)
            && is_finite(qpi_input0)
            && is_finite(xqpi0)
            && is_finite(q_aw_rate0)
            && is_finite(xqpi_rate0)
            && is_finite(ef0)
            && is_finite(pfreq0)
            && is_finite(ep0)
            && is_finite(eplim0)
            && is_finite(pref0)
            && is_finite(ppi0)
            && is_finite(ppi_input0)
            && is_finite(xppi0)
            && is_finite(p_aw_rate0)
            && is_finite(pext_output0)
            && is_finite(freqref0)
            && is_finite(vref0)
            && is_finite(qref0_system)
            && is_finite(pref0_system);
        if (!candidates_are_finite)
        {
          Log::error() << "Repca: derived initial values must be finite\n";
          return 1;
        }

        y[VMEAS]  = vmeas0;
        y[QMEAS]  = qmeas0;
        y[XQPI]   = xqpi0;
        y[XQLAG]  = xqlag0;
        y[PMEAS]  = pmeas0;
        y[XPPI]   = xppi0;
        y[PREF]   = pref0;
        y[V]      = v0;
        y[VLDC]   = vldc0;
        y[VDROOP] = vdroop0;
        y[VCTRL]  = vctrl0;
        y[SFRZ]   = sfrz0;
        y[ERQ]    = erq0;
        y[ERQDB]  = erqdb0;
        y[ERQLIM] = erqlim0;
        y[QPI]    = qpi0;
        y[QEXT]   = qext0_system;
        y[EF]     = ef0;
        y[EP]     = ep0;
        y[EPLIM]  = eplim0;
        y[PPI]    = ppi0;
        y[PEXT]   = pext_output0;

        const bool q_adjusted = qmin != Qmin_ || qmax != Qmax_;
        const bool p_adjusted = pmin != Pmin_ || pmax != Pmax_;
        Qmin_                 = qmin;
        Qmax_                 = qmax;
        Pmin_                 = pmin;
        Pmax_                 = pmax;

        freqref_set_ = freqref0;
        vref_set_    = vref0;
        qref_set_    = qref0_system;
        pref_set_    = pref0_system;

        if (signals_.template isAttached<RepcaExternalVariables::VREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::VREF>(vref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::PREF>(
              pref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::QREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::QREF>(qref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::FREQREF>(
              freqref_set_);
        }

        if (q_adjusted)
        {
          Log::warning() << "Repca: initial reactive PI output is outside [Qmin, Qmax]; "
                            "the limits are adjusted to include the initialized value\n";
        }
        if (p_adjusted)
        {
          Log::warning() << "Repca: initial active PI output is outside [Pmin, Pmax]; "
                            "the limits are adjusted to include the initialized value\n";
        }
        if (Freqflag_
            && (Ddn_ != ZERO<RealT> || Dup_ != ZERO<RealT>)
            && !signals_.template isAttached<RepcaExternalVariables::FREQ>())
        {
          Log::warning() << "Repca: Freqflag is enabled without a freq signal; "
                            "frequency remains nominal and cannot track bus-frequency deviations\n";
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief Identify the differential variables
       *
       * The voltage and power filters, reactive PI and lead-lag states, and
       * active PI and command-lag states carry derivatives; every other
       * internal variable is algebraic.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::tagDifferentiable()
      {
        const auto VMEAS = static_cast<size_t>(RepcaInternalVariables::VMEAS);
        const auto QMEAS = static_cast<size_t>(RepcaInternalVariables::QMEAS);
        const auto XQPI  = static_cast<size_t>(RepcaInternalVariables::XQPI);
        const auto XQLAG = static_cast<size_t>(RepcaInternalVariables::XQLAG);
        const auto PMEAS = static_cast<size_t>(RepcaInternalVariables::PMEAS);
        const auto XPPI  = static_cast<size_t>(RepcaInternalVariables::XPPI);
        const auto PREF  = static_cast<size_t>(RepcaInternalVariables::PREF);

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[VMEAS] = true;
        tag_[QMEAS] = true;
        tag_[XQPI]  = true;
        tag_[XQLAG] = true;
        tag_[PMEAS] = true;
        tag_[XPPI]  = true;
        tag_[PREF]  = true;
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * REPCA variables are dimensionless per-unit signals and controller
       * states of the same order, so they share the relative tolerance as
       * their absolute floor.
       *
       * @param[in] rel_tol Solver relative tolerance.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Evaluate the REPCA-owned residual rows
       *
       * Refreshes required bus and measurement inputs, starts optional
       * references from the values latched by initialize(), then overwrites
       * them from attached signals. REPCA contributes no bus residual.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateResidual()
      {
        const auto IR      = static_cast<size_t>(RepcaExternalVariables::IR);
        const auto II      = static_cast<size_t>(RepcaExternalVariables::II);
        const auto P       = static_cast<size_t>(RepcaExternalVariables::P);
        const auto Q       = static_cast<size_t>(RepcaExternalVariables::Q);
        const auto FREQ    = static_cast<size_t>(RepcaExternalVariables::FREQ);
        const auto VREF    = static_cast<size_t>(RepcaExternalVariables::VREF);
        const auto PREF    = static_cast<size_t>(RepcaExternalVariables::PREF);
        const auto QREF    = static_cast<size_t>(RepcaExternalVariables::QREF);
        const auto FREQREF = static_cast<size_t>(RepcaExternalVariables::FREQREF);

        ws_[VREF]    = vref_set_;
        ws_[PREF]    = pref_set_;
        ws_[QREF]    = qref_set_;
        ws_[FREQ]    = static_cast<ScalarT>(ONE<RealT>);
        ws_[FREQREF] = freqref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        ws_[IR] =
            signals_.template readExternalVariable<RepcaExternalVariables::IR>();
        ws_indices_[IR] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::IR>();
        ws_[II] =
            signals_.template readExternalVariable<RepcaExternalVariables::II>();
        ws_indices_[II] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::II>();
        ws_[P] =
            signals_.template readExternalVariable<RepcaExternalVariables::P>();
        ws_indices_[P] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::P>();
        ws_[Q] =
            signals_.template readExternalVariable<RepcaExternalVariables::Q>();
        ws_indices_[Q] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::Q>();
        if (signals_.template isAttached<RepcaExternalVariables::FREQ>())
        {
          ws_[FREQ] =
              signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();
          ws_indices_[FREQ] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQ>();
        }

        if (signals_.template isAttached<RepcaExternalVariables::VREF>())
        {
          ws_[VREF] =
              signals_.template readExternalVariable<RepcaExternalVariables::VREF>();
          ws_indices_[VREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::VREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::PREF>())
        {
          ws_[PREF] =
              signals_.template readExternalVariable<RepcaExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::PREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::QREF>())
        {
          ws_[QREF] =
              signals_.template readExternalVariable<RepcaExternalVariables::QREF>();
          ws_indices_[QREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::QREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
        {
          ws_[FREQREF] =
              signals_.template readExternalVariable<RepcaExternalVariables::FREQREF>();
          ws_indices_[FREQREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQREF>();
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
       * @brief Access the REPCA signal interface
       *
       * @return Interface used to assign optional plant-command outputs and
       *         attach required measurements and optional references.
       */
      template <typename scalar_type, typename index_type>
      auto Repca<scalar_type, index_type>::getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              RepcaInternalVariables,
                              RepcaExternalVariables>&
      {
        return signals_;
      }

      /**
       * @brief Access the configured monitor
       *
       * @return Monitor for this model, or nullptr when constructed without data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Repca<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Evaluate the REPCA internal residual
       *
       * Evaluates seven differential and fifteen algebraic rows in enum order.
       * Precomputed mode masks keep the differentiated path branch- and
       * loop-free with a fixed dependency structure.
       *
       * @param[in] y Internal variables in `RepcaInternalVariables` order and
       *              on the bases documented by their enums.
       * @param[in] yp Internal derivatives in the same enum order and bases.
       * @param[in] wb Regulated-bus real and imaginary voltage components.
       * @param[in] ws External signals in `RepcaExternalVariables` order.
       * @param[out] f Caller-provided residual output buffer in
       *               `RepcaInternalVariables` order.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline int
      Repca<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto VMEAS      = static_cast<size_t>(RepcaInternalVariables::VMEAS);
        const auto QMEAS      = static_cast<size_t>(RepcaInternalVariables::QMEAS);
        const auto XQPI       = static_cast<size_t>(RepcaInternalVariables::XQPI);
        const auto XQLAG      = static_cast<size_t>(RepcaInternalVariables::XQLAG);
        const auto PMEAS      = static_cast<size_t>(RepcaInternalVariables::PMEAS);
        const auto XPPI       = static_cast<size_t>(RepcaInternalVariables::XPPI);
        const auto PREF_STATE = static_cast<size_t>(RepcaInternalVariables::PREF);
        const auto V          = static_cast<size_t>(RepcaInternalVariables::V);
        const auto VLDC       = static_cast<size_t>(RepcaInternalVariables::VLDC);
        const auto VDROOP     = static_cast<size_t>(RepcaInternalVariables::VDROOP);
        const auto VCTRL      = static_cast<size_t>(RepcaInternalVariables::VCTRL);
        const auto SFRZ       = static_cast<size_t>(RepcaInternalVariables::SFRZ);
        const auto ERQ        = static_cast<size_t>(RepcaInternalVariables::ERQ);
        const auto ERQDB      = static_cast<size_t>(RepcaInternalVariables::ERQDB);
        const auto ERQLIM     = static_cast<size_t>(RepcaInternalVariables::ERQLIM);
        const auto QPI        = static_cast<size_t>(RepcaInternalVariables::QPI);
        const auto QEXT       = static_cast<size_t>(RepcaInternalVariables::QEXT);
        const auto EF         = static_cast<size_t>(RepcaInternalVariables::EF);
        const auto EP         = static_cast<size_t>(RepcaInternalVariables::EP);
        const auto EPLIM      = static_cast<size_t>(RepcaInternalVariables::EPLIM);
        const auto PPI        = static_cast<size_t>(RepcaInternalVariables::PPI);
        const auto PEXT       = static_cast<size_t>(RepcaInternalVariables::PEXT);

        const auto IR         = static_cast<size_t>(RepcaExternalVariables::IR);
        const auto II         = static_cast<size_t>(RepcaExternalVariables::II);
        const auto P          = static_cast<size_t>(RepcaExternalVariables::P);
        const auto Q          = static_cast<size_t>(RepcaExternalVariables::Q);
        const auto FREQ       = static_cast<size_t>(RepcaExternalVariables::FREQ);
        const auto VREF       = static_cast<size_t>(RepcaExternalVariables::VREF);
        const auto PREF_INPUT = static_cast<size_t>(RepcaExternalVariables::PREF);
        const auto QREF       = static_cast<size_t>(RepcaExternalVariables::QREF);
        const auto FREQREF    = static_cast<size_t>(RepcaExternalVariables::FREQREF);

        const ScalarT vmeas  = y[VMEAS];
        const ScalarT qmeas  = y[QMEAS];
        const ScalarT xqpi   = y[XQPI];
        const ScalarT xqlag  = y[XQLAG];
        const ScalarT pmeas  = y[PMEAS];
        const ScalarT xppi   = y[XPPI];
        const ScalarT pref   = y[PREF_STATE];
        const ScalarT v      = y[V];
        const ScalarT vldc   = y[VLDC];
        const ScalarT vdroop = y[VDROOP];
        const ScalarT vctrl  = y[VCTRL];
        const ScalarT sfrz   = y[SFRZ];
        const ScalarT erq    = y[ERQ];
        const ScalarT erqdb  = y[ERQDB];
        const ScalarT erqlim = y[ERQLIM];
        const ScalarT qpi    = y[QPI];
        const ScalarT qext   = toComponentBase(y[QEXT]);
        const ScalarT ef     = y[EF];
        const ScalarT ep     = y[EP];
        const ScalarT eplim  = y[EPLIM];
        const ScalarT ppi    = y[PPI];
        const ScalarT pext   = toComponentBase(y[PEXT]);

        const ScalarT vmeas_dot = yp[VMEAS];
        const ScalarT qmeas_dot = yp[QMEAS];
        const ScalarT xqpi_dot  = yp[XQPI];
        const ScalarT xqlag_dot = yp[XQLAG];
        const ScalarT pmeas_dot = yp[PMEAS];
        const ScalarT xppi_dot  = yp[XPPI];
        const ScalarT pref_dot  = yp[PREF_STATE];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT ir      = toComponentBase(ws[IR]);
        const ScalarT ii      = toComponentBase(ws[II]);
        const ScalarT p       = toComponentBase(ws[P]);
        const ScalarT q       = toComponentBase(ws[Q]);
        const ScalarT freq    = ws[FREQ];
        const ScalarT freqref = ws[FREQREF];
        const ScalarT vref    = ws[VREF];
        const ScalarT qref    = toComponentBase(ws[QREF]);
        const ScalarT pref_in = toComponentBase(ws[PREF_INPUT]);

        const ScalarT vldc_r = vr - Rc_ * ir + Xc_ * ii;
        const ScalarT vldc_i = vi - Rc_ * ii - Xc_ * ir;
        const ScalarT pfreq  = droop(ef, Ddn_, Dup_);

        f[VMEAS]      = -vmeas_dot + (vctrl - vmeas) / Tfltr_;
        f[QMEAS]      = -qmeas_dot + (q - qmeas) / Tfltr_;
        f[XQPI]       = -xqpi_dot + sfrz * Math::antiwindup(qpi, Ki_ * erqlim, Qmin_, Qmax_);
        f[XQLAG]      = -xqlag_dot + (qpi - xqlag) / Tfv_;
        f[PMEAS]      = -pmeas_dot + (p - pmeas) / Tp_;
        f[XPPI]       = -xppi_dot + Math::antiwindup(ppi, Kig_ * eplim, Pmin_, Pmax_);
        f[PREF_STATE] = -pref_dot + (ppi - pref) / Tlag_;

        f[V]      = -v * v + vr * vr + vi * vi;
        f[VLDC]   = -vldc * vldc + vldc_r * vldc_r + vldc_i * vldc_i;
        f[VDROOP] = -vdroop + v + Kc_ * q;
        f[VCTRL]  = -vctrl + vcomp_on_ * vldc + vcomp_off_ * vdroop;
        f[SFRZ]   = -sfrz + Math::above(v, Vfrz_);
        f[ERQ]    = -erq + ref_on_ * (vref - vmeas) + ref_off_ * (qref - qmeas);
        f[ERQDB]  = -erqdb + Math::deadband2(erq, dbdlow_, dbdupper_);
        f[ERQLIM] = -erqlim + Math::clamp(erqdb, emin_, emax_);
        f[QPI]    = -qpi + Math::clamp(Kp_ * erqlim + xqpi, Qmin_, Qmax_);
        f[QEXT]   = -Tfv_ * (qext - xqlag) + Tft_ * (qpi - xqlag);

        f[EF]    = -ef + Math::deadband2(freqref - freq, fdbd1_, fdbd2_);
        f[EP]    = -ep + pref_in - pmeas + pfreq;
        f[EPLIM] = -eplim + Math::clamp(ep, femin_, femax_);
        f[PPI]   = -ppi + Math::clamp(Kpg_ * eplim + xppi, Pmin_, Pmax_);
        f[PEXT]  = -pext + freq_on_ * pref;

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Smooth asymmetric frequency-droop response
       *
       * @param[in] error Deadbanded frequency error.
       * @param[in] down Down-regulation (overfrequency) gain.
       * @param[in] up Up-regulation (underfrequency) gain.
       * @return Active-power frequency response.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline scalar_type
      Repca<scalar_type, index_type>::droop(ScalarT error, RealT down, RealT up)
      {
        return error * (down + (up - down) * Math::sigmoid(error));
      }

      /**
       * @brief Read optional parameters from model data
       *
       * Omitted parameters retain their documented defaults. Real parameters
       * accept real and integer values and must be finite; selectors require
       * Boolean values. Loading errors are counted for verify() rather than
       * thrown.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
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
          RealT       parsed_value{};
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            parsed_value = *real_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            parsed_value = static_cast<RealT>(*index_value);
          }
          else
          {
            Log::error() << "Repca: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
            return;
          }

          if (!std::isfinite(parsed_value))
          {
            Log::error() << "Repca: parameter '" << name << "' must be finite\n";
            ++parameter_error_count_;
            return;
          }

          target = parsed_value;
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
          else
          {
            Log::error() << "Repca: parameter '" << name << "' must be boolean\n";
            ++parameter_error_count_;
          }
        };

        load_real(Params::mva, mva_base_, "mva");
        load_switch(Params::VcompFlag, VcompFlag_, "VcompFlag");
        load_switch(Params::RefFlag, RefFlag_, "RefFlag");
        load_switch(Params::Freqflag, Freqflag_, "Freqflag");
        load_real(Params::Tfltr, Tfltr_, "Tfltr");
        load_real(Params::Vfrz, Vfrz_, "Vfrz");
        load_real(Params::Rc, Rc_, "Rc");
        load_real(Params::Xc, Xc_, "Xc");
        load_real(Params::Kc, Kc_, "Kc");
        load_real(Params::dbdlow, dbdlow_, "dbdlow");
        load_real(Params::dbdupper, dbdupper_, "dbdupper");
        load_real(Params::emax, emax_, "emax");
        load_real(Params::emin, emin_, "emin");
        load_real(Params::Kp, Kp_, "Kp");
        load_real(Params::Ki, Ki_, "Ki");
        load_real(Params::Qmax, Qmax_, "Qmax");
        load_real(Params::Qmin, Qmin_, "Qmin");
        load_real(Params::Tft, Tft_, "Tft");
        load_real(Params::Tfv, Tfv_, "Tfv");
        load_real(Params::Tp, Tp_, "Tp");
        load_real(Params::fdbd1, fdbd1_, "fdbd1");
        load_real(Params::fdbd2, fdbd2_, "fdbd2");
        load_real(Params::Ddn, Ddn_, "Ddn");
        load_real(Params::Dup, Dup_, "Dup");
        load_real(Params::femax, femax_, "femax");
        load_real(Params::femin, femin_, "femin");
        load_real(Params::Kpg, Kpg_, "Kpg");
        load_real(Params::Kig, Kig_, "Kig");
        load_real(Params::Pmax, Pmax_, "Pmax");
        load_real(Params::Pmin, Pmin_, "Pmin");
        load_real(Params::Tlag, Tlag_, "Tlag");

        setDerivedParameters();
      }

      /**
       * @brief Bind monitor selections to REPCA internal states
       */
      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::qext, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::QEXT)]; });
        monitor_->set(Variable::pext, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::PEXT)]; });
        monitor_->set(Variable::vmeas, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::VMEAS)]; });
        monitor_->set(Variable::qmeas, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::QMEAS)]; });
        monitor_->set(Variable::pmeas, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::PMEAS)]; });
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "Repca: Tfltr, Tfv, Tp, and Tlag below "
                       << TIME_CONSTANT_MINIMUM
                       << " s are raised to that floor to keep the controller lags well posed\n";
      }

      /**
       * @brief Resolve parameter-derived constants
       *
       * Raises the explicit controller lags in place, computes the component
       * power base, and resolves selector masks.
       */
      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::setDerivedParameters()
      {
        // The lags are raised to the floor below, so negative values must be
        // rejected here while the value as read is still available.
        auto check_non_negative = [&](RealT value, const char* name)
        {
          if (value < ZERO<RealT>)
          {
            Log::error() << "Repca: " << name << " must be non-negative\n";
            ++parameter_error_count_;
          }
        };

        check_non_negative(Tfltr_, "Tfltr");
        check_non_negative(Tft_, "Tft");
        check_non_negative(Tfv_, "Tfv");
        check_non_negative(Tp_, "Tp");
        check_non_negative(Tlag_, "Tlag");

        if (Tfltr_ < TIME_CONSTANT_MINIMUM || Tfv_ < TIME_CONSTANT_MINIMUM
            || Tp_ < TIME_CONSTANT_MINIMUM || Tlag_ < TIME_CONSTANT_MINIMUM)
        {
          static std::once_flag time_constant_warning_flag_;
          std::call_once(time_constant_warning_flag_,
                         &logTimeConstantWarning);
        }

        Tfltr_ = std::max(Tfltr_, TIME_CONSTANT_MINIMUM);
        Tfv_   = std::max(Tfv_, TIME_CONSTANT_MINIMUM);
        Tp_    = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        Tlag_  = std::max(Tlag_, TIME_CONSTANT_MINIMUM);

        va_component_base_ = mva_base_ * static_cast<RealT>(1.0e6);

        vcomp_on_  = VcompFlag_ ? ONE<RealT> : ZERO<RealT>;
        vcomp_off_ = ONE<RealT> - vcomp_on_;
        ref_on_    = RefFlag_ ? ONE<RealT> : ZERO<RealT>;
        ref_off_   = ONE<RealT> - ref_on_;
        freq_on_   = Freqflag_ ? ONE<RealT> : ZERO<RealT>;
      }

      /**
       * @brief Recover an input for a requested smooth-clamp output
       *
       * This initialization-only helper inverts the CommonMath smooth clamp,
       * including collapsed limits. Exact-bound requests use a 0.1 offset,
       * leaving less than 2e-13 at CommonMath MU = 240.
       *
       * @param[in] output Requested smooth-clamp output.
       * @param[in] lower Lower clamp limit.
       * @param[in] upper Upper clamp limit.
       * @param[out] input Recovered clamp input on success.
       * @return true when a finite admissible input was recovered.
       * @warning Contains conditional branching and is for initialization only;
       *          do not use it during residual evaluation.
       */
      template <typename scalar_type, typename index_type>
      bool Repca<scalar_type, index_type>::invertClamp(ScalarT output, RealT lower, RealT upper, ScalarT& input) const
      {
        const RealT value = static_cast<RealT>(output);

        if (!std::isfinite(value)
            || !std::isfinite(lower)
            || !std::isfinite(upper)
            || lower > upper
            || value < lower
            || value > upper)
        {
          return false;
        }

        const RealT width = upper - lower;
        if (width <= INITIALIZATION_TOLERANCE)
        {
          input = static_cast<ScalarT>(lower);
          return true;
        }

        const RealT distance_from_lower = value - lower;
        const RealT distance_from_upper = upper - value;
        if (distance_from_lower <= INITIALIZATION_TOLERANCE)
        {
          input = static_cast<ScalarT>(lower - INITIALIZATION_LIMIT_OFFSET);
          return true;
        }
        if (distance_from_upper <= INITIALIZATION_TOLERANCE)
        {
          input = static_cast<ScalarT>(upper + INITIALIZATION_LIMIT_OFFSET);
          return true;
        }

        const RealT mu                    = Math::MU<RealT>;
        const RealT scaled_lower_distance = mu * distance_from_lower;
        const RealT scaled_upper_distance = mu * distance_from_upper;
        const RealT log_lower             = logOneMinusExp(scaled_lower_distance);
        const RealT log_upper             = logOneMinusExp(scaled_upper_distance);
        const RealT correction            = (scaled_lower_distance + log_lower - log_upper) / mu;

        input = static_cast<ScalarT>(lower + correction);
        return std::isfinite(static_cast<RealT>(input));
      }

      /**
       * @brief Recover an input for a requested smooth deadband output
       *
       * @param[in] output Requested deadband output.
       * @param[in] lower Lower deadband threshold.
       * @param[in] upper Upper deadband threshold.
       * @param[out] input Recovered deadband input on success.
       * @return true when a finite input was recovered.
       * @warning Contains conditional branching and is for initialization only;
       *          do not use it during residual evaluation.
       */
      template <typename scalar_type, typename index_type>
      bool Repca<scalar_type, index_type>::invertDeadband(ScalarT output, RealT lower, RealT upper, ScalarT& input) const
      {
        const RealT value = static_cast<RealT>(output);

        if (!std::isfinite(value)
            || !std::isfinite(lower)
            || !std::isfinite(upper)
            || lower > upper)
        {
          return false;
        }

        const RealT midpoint = HALF<RealT> * lower + HALF<RealT> * upper;
        if (std::abs(value) <= INITIALIZATION_TOLERANCE)
        {
          input = static_cast<ScalarT>(midpoint);
          return true;
        }

        RealT lower_input = midpoint;
        RealT upper_input = midpoint;
        if (value < ZERO<RealT>)
        {
          lower_input = lower + value;
        }
        else
        {
          upper_input = upper + value;
        }

        const RealT lower_output = Math::deadband2(lower_input, lower, upper);
        const RealT upper_output = Math::deadband2(upper_input, lower, upper);
        if (!std::isfinite(lower_output)
            || !std::isfinite(upper_output)
            || lower_output - value > INITIALIZATION_TOLERANCE
            || value - upper_output > INITIALIZATION_TOLERANCE)
        {
          return false;
        }

        for (std::size_t iteration = 0; iteration < 128; ++iteration)
        {
          const RealT mid = lower_input + HALF<RealT> * (upper_input - lower_input);
          if (Math::deadband2(mid, lower, upper) < value)
          {
            lower_input = mid;
          }
          else
          {
            upper_input = mid;
          }
        }

        const RealT result = lower_input + HALF<RealT> * (upper_input - lower_input);
        input              = static_cast<ScalarT>(result);
        return std::isfinite(result)
               && std::abs(Math::deadband2(result, lower, upper) - value)
                      <= INITIALIZATION_TOLERANCE;
      }

      /**
       * @brief Evaluate log(1 - exp(-x)) accurately for positive x
       *
       * The small-x form avoids cancellation in `1 - exp(-x)`; the large-x
       * form uses `log1p`. The two algebraically equivalent forms agree in
       * value and first derivative at x = log(2).
       *
       * @param[in] x Positive argument.
       * @return Numerically stable value of log(1 - exp(-x)).
       */
      template <typename scalar_type, typename index_type>
      typename Repca<scalar_type, index_type>::RealT
      Repca<scalar_type, index_type>::logOneMinusExp(RealT x)
      {
        static constexpr auto log_two = std::numbers::ln2_v<RealT>;

        if (x < log_two)
        {
          return log_two - HALF<RealT> * x
                 + std::log(std::sinh(HALF<RealT> * x));
        }
        return std::log1p(-std::exp(-x));
      }

      /**
       * @brief Convert a system-base power or current quantity to REPCA component base
       *
       * @param[in] value Quantity on the system base.
       * @return The same quantity on the REPCA component base.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      Repca<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * (va_system_base_ / va_component_base_);
      }

      /**
       * @brief Convert a component-base power quantity to system base
       *
       * @param[in] value Quantity on the REPCA component base.
       * @return The same quantity on the system base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Repca<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value * (va_component_base_ / va_system_base_);
      }

      /**
       * @brief Access the regulated-bus real voltage component
       *
       * @return Mutable reference to the bus real voltage state.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Repca<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      /**
       * @brief Access the regulated-bus imaginary voltage component
       *
       * @return Mutable reference to the bus imaginary voltage state.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Repca<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit

/**
 * @file RepcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REPCA phasor-dynamics plant-control model.
 */

#pragma once

#include <algorithm>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Logger used for REPCA diagnostics.
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @class Repca
       * @brief WECC renewable plant controller with reactive-power and
       *        active-power control paths.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */

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
        size_ = static_cast<IdxT>(index(RepcaInternalVariables::MAXIMUM));
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
        size_ = static_cast<IdxT>(index(RepcaInternalVariables::MAXIMUM));
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
        using I = RepcaInternalVariables;
        using E = RepcaExternalVariables;

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        const auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(2, ScalarT{0});

        const auto signal_size = index(E::MAXIMUM);
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
              &y[index(I::QEXT)],
              &(this->getVariableIndex(static_cast<IdxT>(index(I::QEXT)))));
        }

        if (signals_.template isAssigned<RepcaInternalVariables::PEXT>())
        {
          signals_.template getSignalNode<RepcaInternalVariables::PEXT>()->set(
              &y[index(I::PEXT)],
              &(this->getVariableIndex(static_cast<IdxT>(index(I::PEXT)))));
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

        const bool valid_converter_base = std::isfinite(mva_base_)
                                          && mva_base_ > ZERO<RealT>
                                          && std::isfinite(va_converter_base_)
                                          && va_converter_base_ > ZERO<RealT>;
        const bool valid_system_base = std::isfinite(va_system_base_)
                                       && va_system_base_ > ZERO<RealT>;
        check(valid_converter_base,
              "mva must define a finite positive converter power base");
        check(valid_system_base, "system power base must be finite and positive");
        if (valid_converter_base && valid_system_base)
        {
          const RealT system_to_converter = va_system_base_ / va_converter_base_;
          const RealT converter_to_system = va_converter_base_ / va_system_base_;
          check(std::isfinite(system_to_converter)
                    && system_to_converter > ZERO<RealT>
                    && std::isfinite(converter_to_system)
                    && converter_to_system > ZERO<RealT>,
                "system/converter power-base conversion ratios must be finite and positive");
        }

        check(dbdlow_ <= ZERO<RealT> && ZERO<RealT> <= dbdupper_,
              "dbdlow <= 0 <= dbdupper is required");
        check(emin_ <= ZERO<RealT> && ZERO<RealT> <= emax_,
              "emin <= 0 <= emax is required");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(Tfv_ > ZERO<RealT>, "Tfv must be positive");
        check(fdbd1_ <= ZERO<RealT> && ZERO<RealT> <= fdbd2_,
              "fdbd1 <= 0 <= fdbd2 is required");
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
        check_required_signal.template operator()<RepcaExternalVariables::FREQ>("freq");

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
       * @post On failure no state, derivative, latch, or signal storage is
       *       modified.
       *
       * @return 0 on success; nonzero when allocation, configuration, initial-
       *         value, candidate, command-limit, or steady-state checks fail.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::initialize()
      {
        using I = RepcaInternalVariables;

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

        const ScalarT qext0_system = y[index(I::QEXT)];
        const ScalarT pext0_system = y[index(I::PEXT)];
        const ScalarT qext0        = toComponentBase(qext0_system);
        const ScalarT pext0        = toComponentBase(pext0_system);

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();
        const ScalarT ir =
            signals_.template readExternalVariable<RepcaExternalVariables::IR>();
        const ScalarT ii =
            signals_.template readExternalVariable<RepcaExternalVariables::II>();
        const ScalarT p_system =
            signals_.template readExternalVariable<RepcaExternalVariables::P>();
        const ScalarT q_system =
            signals_.template readExternalVariable<RepcaExternalVariables::Q>();
        const ScalarT freq = signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();

        auto is_finite = [](ScalarT value)
        {
          return std::isfinite(static_cast<RealT>(value));
        };
        if (!is_finite(vr) || !is_finite(vi) || !is_finite(ir)
            || !is_finite(ii) || !is_finite(p_system)
            || !is_finite(q_system)
            || !is_finite(freq) || !is_finite(qext0)
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

        const ScalarT zero    = static_cast<ScalarT>(ZERO<RealT>);
        const ScalarT erq0    = zero;
        const ScalarT erqdb0  = Math::deadband2(erq0, dbdlow_, dbdupper_);
        const ScalarT erqlim0 = Math::clamp(erqdb0, emin_, emax_);
        const ScalarT qpi0    = qext0;
        const ScalarT xqlag0  = qpi0;

        ScalarT qpi_input0{};
        if (!solveLimiterInput(qpi0, Qmin_, Qmax_, qpi_input0))
        {
          Log::error() << "Repca: initial reactive-power command is outside Qmin/Qmax\n";
          return 1;
        }
        const ScalarT xqpi0      = qpi_input0 - Kp_ * erqlim0;
        const ScalarT q_aw_rate0 = Math::antiwindup(qpi0, Ki_ * erqlim0, Qmin_, Qmax_);
        const ScalarT xqpi_rate0 = sfrz0 * q_aw_rate0;
        if (!is_finite(q_aw_rate0) || !is_finite(xqpi_rate0)
            || std::abs(static_cast<RealT>(xqpi_rate0)) > INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Repca: reactive-power PI state rate is nonzero at initialization\n";
          return 1;
        }

        const ScalarT ef0    = Math::deadband2(zero, fdbd1_, fdbd2_);
        const ScalarT pfreq0 = Ddn_ * Math::ramp(ef0) - Dup_ * Math::ramp(-ef0);
        const ScalarT ep0    = zero;
        const ScalarT eplim0 = Math::clamp(ep0, femin_, femax_);
        const ScalarT pref0  = Freqflag_ ? pext0 : Math::clamp(pmeas0, Pmin_, Pmax_);
        const ScalarT ppi0   = pref0;
        ScalarT       ppi_input0{};
        if (Freqflag_)
        {
          if (!solveLimiterInput(ppi0, Pmin_, Pmax_, ppi_input0))
          {
            Log::error() << "Repca: initial active-power command is outside Pmin/Pmax\n";
            return 1;
          }
        }
        else
        {
          ppi_input0 = pmeas0;
        }
        const ScalarT xppi0      = ppi_input0 - Kpg_ * eplim0;
        const ScalarT p_aw_rate0 = Math::antiwindup(ppi0, Kig_ * eplim0, Pmin_, Pmax_);
        if (!is_finite(p_aw_rate0)
            || std::abs(static_cast<RealT>(p_aw_rate0)) > INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Repca: active-power PI antiwindup rate is nonzero at initialization\n";
          return 1;
        }

        const ScalarT pext_output0 = Freqflag_ ? pext0_system : zero;
        const ScalarT freqref0     = freq;
        const ScalarT vref0        = vmeas0;
        const ScalarT qref0_system = q_system;
        const ScalarT pref0_system = toSystemBase(pmeas0 - pfreq0);

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

        y[index(I::VMEAS)]  = vmeas0;
        y[index(I::QMEAS)]  = qmeas0;
        y[index(I::XQPI)]   = xqpi0;
        y[index(I::XQLAG)]  = xqlag0;
        y[index(I::PMEAS)]  = pmeas0;
        y[index(I::XPPI)]   = xppi0;
        y[index(I::PREF)]   = pref0;
        y[index(I::V)]      = v0;
        y[index(I::VLDC)]   = vldc0;
        y[index(I::VDROOP)] = vdroop0;
        y[index(I::VCTRL)]  = vctrl0;
        y[index(I::SFRZ)]   = sfrz0;
        y[index(I::ERQ)]    = erq0;
        y[index(I::ERQDB)]  = erqdb0;
        y[index(I::ERQLIM)] = erqlim0;
        y[index(I::QPI)]    = qpi0;
        y[index(I::QEXT)]   = qext0_system;
        y[index(I::EF)]     = ef0;
        y[index(I::EP)]     = ep0;
        y[index(I::EPLIM)]  = eplim0;
        y[index(I::PPI)]    = ppi0;
        y[index(I::PEXT)]   = pext_output0;

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
        using I = RepcaInternalVariables;

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[index(I::VMEAS)] = true;
        tag_[index(I::QMEAS)] = true;
        tag_[index(I::XQPI)]  = true;
        tag_[index(I::XQLAG)] = true;
        tag_[index(I::PMEAS)] = true;
        tag_[index(I::XPPI)]  = true;
        tag_[index(I::PREF)]  = true;
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
        using E = RepcaExternalVariables;

        ws_[index(E::VREF)]    = vref_set_;
        ws_[index(E::PREF)]    = pref_set_;
        ws_[index(E::QREF)]    = qref_set_;
        ws_[index(E::FREQREF)] = freqref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        ws_[index(E::IR)] =
            signals_.template readExternalVariable<RepcaExternalVariables::IR>();
        ws_indices_[index(E::IR)] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::IR>();
        ws_[index(E::II)] =
            signals_.template readExternalVariable<RepcaExternalVariables::II>();
        ws_indices_[index(E::II)] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::II>();
        ws_[index(E::P)] =
            signals_.template readExternalVariable<RepcaExternalVariables::P>();
        ws_indices_[index(E::P)] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::P>();
        ws_[index(E::Q)] =
            signals_.template readExternalVariable<RepcaExternalVariables::Q>();
        ws_indices_[index(E::Q)] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::Q>();
        ws_[index(E::FREQ)] =
            signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();
        ws_indices_[index(E::FREQ)] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQ>();

        if (signals_.template isAttached<RepcaExternalVariables::VREF>())
        {
          ws_[index(E::VREF)] =
              signals_.template readExternalVariable<RepcaExternalVariables::VREF>();
          ws_indices_[index(E::VREF)] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::VREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::PREF>())
        {
          ws_[index(E::PREF)] =
              signals_.template readExternalVariable<RepcaExternalVariables::PREF>();
          ws_indices_[index(E::PREF)] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::PREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::QREF>())
        {
          ws_[index(E::QREF)] =
              signals_.template readExternalVariable<RepcaExternalVariables::QREF>();
          ws_indices_[index(E::QREF)] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::QREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
        {
          ws_[index(E::FREQREF)] =
              signals_.template readExternalVariable<RepcaExternalVariables::FREQREF>();
          ws_indices_[index(E::FREQREF)] =
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
       * @param[out] f Model-owned residuals in `RepcaInternalVariables` order.
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
        using I = RepcaInternalVariables;
        using E = RepcaExternalVariables;

        const ScalarT vmeas  = y[index(I::VMEAS)];
        const ScalarT qmeas  = y[index(I::QMEAS)];
        const ScalarT xqpi   = y[index(I::XQPI)];
        const ScalarT xqlag  = y[index(I::XQLAG)];
        const ScalarT pmeas  = y[index(I::PMEAS)];
        const ScalarT xppi   = y[index(I::XPPI)];
        const ScalarT pref   = y[index(I::PREF)];
        const ScalarT v      = y[index(I::V)];
        const ScalarT vldc   = y[index(I::VLDC)];
        const ScalarT vdroop = y[index(I::VDROOP)];
        const ScalarT vctrl  = y[index(I::VCTRL)];
        const ScalarT sfrz   = y[index(I::SFRZ)];
        const ScalarT erq    = y[index(I::ERQ)];
        const ScalarT erqdb  = y[index(I::ERQDB)];
        const ScalarT erqlim = y[index(I::ERQLIM)];
        const ScalarT qpi    = y[index(I::QPI)];
        const ScalarT qext   = toComponentBase(y[index(I::QEXT)]);
        const ScalarT ef     = y[index(I::EF)];
        const ScalarT ep     = y[index(I::EP)];
        const ScalarT eplim  = y[index(I::EPLIM)];
        const ScalarT ppi    = y[index(I::PPI)];
        const ScalarT pext   = toComponentBase(y[index(I::PEXT)]);

        const ScalarT vmeas_dot = yp[index(I::VMEAS)];
        const ScalarT qmeas_dot = yp[index(I::QMEAS)];
        const ScalarT xqpi_dot  = yp[index(I::XQPI)];
        const ScalarT xqlag_dot = yp[index(I::XQLAG)];
        const ScalarT pmeas_dot = yp[index(I::PMEAS)];
        const ScalarT xppi_dot  = yp[index(I::XPPI)];
        const ScalarT pref_dot  = yp[index(I::PREF)];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT ir      = ws[index(E::IR)];
        const ScalarT ii      = ws[index(E::II)];
        const ScalarT p       = toComponentBase(ws[index(E::P)]);
        const ScalarT q       = toComponentBase(ws[index(E::Q)]);
        const ScalarT freq    = ws[index(E::FREQ)];
        const ScalarT freqref = ws[index(E::FREQREF)];
        const ScalarT vref    = ws[index(E::VREF)];
        const ScalarT qref    = toComponentBase(ws[index(E::QREF)]);
        const ScalarT pref_in = toComponentBase(ws[index(E::PREF)]);

        const ScalarT vldc_r = vr - Rc_ * ir + Xc_ * ii;
        const ScalarT vldc_i = vi - Rc_ * ii - Xc_ * ir;
        const ScalarT pfreq  = Ddn_ * Math::ramp(ef) - Dup_ * Math::ramp(-ef);

        f[index(I::VMEAS)] = -vmeas_dot + (vctrl - vmeas) / Tfltr_;
        f[index(I::QMEAS)] = -qmeas_dot + (q - qmeas) / Tfltr_;
        f[index(I::XQPI)]  = -xqpi_dot + sfrz * Math::antiwindup(qpi, Ki_ * erqlim, Qmin_, Qmax_);
        f[index(I::XQLAG)] = -xqlag_dot + (qpi - xqlag) / Tfv_;
        f[index(I::PMEAS)] = -pmeas_dot + (p - pmeas) / Tp_;
        f[index(I::XPPI)]  = -xppi_dot + Math::antiwindup(ppi, Kig_ * eplim, Pmin_, Pmax_);
        f[index(I::PREF)]  = -pref_dot + (ppi - pref) / Tlag_;

        f[index(I::V)]      = -v * v + vr * vr + vi * vi;
        f[index(I::VLDC)]   = -vldc * vldc + vldc_r * vldc_r + vldc_i * vldc_i;
        f[index(I::VDROOP)] = -vdroop + v + Kc_ * q;
        f[index(I::VCTRL)]  = -vctrl + vcomp_on_ * vldc + vcomp_off_ * vdroop;
        f[index(I::SFRZ)]   = -sfrz + Math::above(v, Vfrz_);
        f[index(I::ERQ)]    = -erq + ref_on_ * (vref - vmeas) + ref_off_ * (qref - qmeas);
        f[index(I::ERQDB)]  = -erqdb + Math::deadband2(erq, dbdlow_, dbdupper_);
        f[index(I::ERQLIM)] = -erqlim + Math::clamp(erqdb, emin_, emax_);
        f[index(I::QPI)]    = -qpi + Math::clamp(Kp_ * erqlim + xqpi, Qmin_, Qmax_);
        f[index(I::QEXT)]   = -Tfv_ * (qext - xqlag) + Tft_ * (qpi - xqlag);

        f[index(I::EF)]    = -ef + Math::deadband2(freqref - freq, fdbd1_, fdbd2_);
        f[index(I::EP)]    = -ep + pref_in - pmeas + pfreq;
        f[index(I::EPLIM)] = -eplim + Math::clamp(ep, femin_, femax_);
        f[index(I::PPI)]   = -ppi + Math::clamp(Kpg_ * eplim + xppi, Pmin_, Pmax_);
        f[index(I::PEXT)]  = -pext + freq_on_ * pref;

        return 0;
      }

      //
      //  Private methods
      //

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
        using I        = RepcaInternalVariables;
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::qext, [this]
                      { return y_.getData()[index(I::QEXT)]; });
        monitor_->set(Variable::pext, [this]
                      { return y_.getData()[index(I::PEXT)]; });
        monitor_->set(Variable::vmeas, [this]
                      { return y_.getData()[index(I::VMEAS)]; });
        monitor_->set(Variable::qmeas, [this]
                      { return y_.getData()[index(I::QMEAS)]; });
        monitor_->set(Variable::pmeas, [this]
                      { return y_.getData()[index(I::PMEAS)]; });
      }

      /**
       * @brief Resolve parameter-derived constants
       *
       * Raises the explicit controller lags in place, computes the converter
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
        check_non_negative(Tp_, "Tp");
        check_non_negative(Tlag_, "Tlag");

        if (Tfltr_ < TIME_CONSTANT_MINIMUM || Tp_ < TIME_CONSTANT_MINIMUM
            || Tlag_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "Repca: Tfltr, Tp, and Tlag below "
                         << TIME_CONSTANT_MINIMUM
                         << " s are raised to that floor to keep the controller lags well posed\n";
        }

        Tfltr_ = std::max(Tfltr_, TIME_CONSTANT_MINIMUM);
        Tp_    = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        Tlag_  = std::max(Tlag_, TIME_CONSTANT_MINIMUM);

        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);

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
       * leaving less than 2e-13 at CommonMath MU = 240. It must not be called
       * from residual evaluation.
       *
       * @param[in] requested_output Requested smooth-clamp output.
       * @param[in] lower_limit Lower clamp limit.
       * @param[in] upper_limit Upper clamp limit.
       * @param[out] limiter_input Recovered clamp input on success.
       * @return true when a finite admissible input was recovered.
       */
      template <typename scalar_type, typename index_type>
      bool Repca<scalar_type, index_type>::solveLimiterInput(
          ScalarT  requested_output,
          RealT    lower_limit,
          RealT    upper_limit,
          ScalarT& limiter_input) const
      {
        const RealT output_value = static_cast<RealT>(requested_output);

        if (!std::isfinite(output_value)
            || !std::isfinite(lower_limit)
            || !std::isfinite(upper_limit)
            || lower_limit > upper_limit
            || output_value < lower_limit
            || output_value > upper_limit)
        {
          return false;
        }

        const RealT width = upper_limit - lower_limit;
        if (width <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(lower_limit);
          return true;
        }

        const RealT distance_from_lower = output_value - lower_limit;
        const RealT distance_from_upper = upper_limit - output_value;
        if (distance_from_lower <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(lower_limit - INITIALIZATION_LIMIT_OFFSET);
          return true;
        }
        if (distance_from_upper <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(upper_limit + INITIALIZATION_LIMIT_OFFSET);
          return true;
        }

        const RealT mu                    = Math::MU<RealT>;
        const RealT scaled_lower_distance = mu * distance_from_lower;
        const RealT scaled_upper_distance = mu * distance_from_upper;
        const RealT log_lower             = logOneMinusExp(scaled_lower_distance);
        const RealT log_upper             = logOneMinusExp(scaled_upper_distance);
        const RealT correction =
            (scaled_lower_distance + log_lower - log_upper) / mu;

        limiter_input = static_cast<ScalarT>(lower_limit + correction);
        return std::isfinite(static_cast<RealT>(limiter_input));
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
      Repca<scalar_type, index_type>::logOneMinusExp(RealT x) const
      {
        const RealT log_two = std::log(static_cast<RealT>(2.0));

        if (x < log_two)
        {
          return log_two - HALF<RealT> * x
                 + std::log(std::sinh(HALF<RealT> * x));
        }
        return std::log1p(-std::exp(-x));
      }

      /**
       * @brief Convert a system-base power quantity to REPCA component base
       *
       * @param[in] value Quantity on the system base.
       * @return The same quantity on the REPCA component base.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      Repca<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * (va_system_base_ / va_converter_base_);
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
        return value * (va_converter_base_ / va_system_base_);
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

    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

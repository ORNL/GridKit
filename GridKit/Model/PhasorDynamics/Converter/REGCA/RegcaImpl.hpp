/**
 * @file RegcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REGCA phasor-dynamics converter model.
 */

#pragma once

#include <algorithm>
#include <mutex>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct a sized but unconfigured REGCA converter.
       *
       * @param[in] bus Terminal bus the converter injects into.
       */
      template <typename scalar_type, typename index_type>
      Regca<scalar_type, index_type>::Regca(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
      }

      /**
       * @brief Construct a REGCA converter from model data.
       *
       * Parameters are loaded and the monitor is created. Model-data errors
       * are recorded for verify() rather than thrown.
       *
       * @param[in] bus Terminal bus the converter injects into.
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      Regca<scalar_type, index_type>::Regca(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Regca<scalar_type, index_type>::~Regca()
      {
      }

      /**
       * @brief Static method to log time constant warnings
       *
       * @note Used in combination with static std:once_flag and std:call_once,
       *       to reduce the number of times the warning is printed.
       */
      template <typename scalar_type, typename index_type>
      void Regca<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "Regca: Tg and TM below " << TIME_CONSTANT_MINIMUM
                       << " s are raised to that floor to keep the current-control "
                       << "and voltage-sensor lags well posed\n";
      }

      /**
       * @brief Resolve parameter-derived constants and limiter selections.
       *
       * The time constants are raised to the well-posedness floor. Complementary
       * LVPL masks keep the residual branchless, while the sign of the initial
       * reactive-power injection selects the applicable recovery-rate limit.
       */
      template <typename scalar_type, typename index_type>
      void Regca<scalar_type, index_type>::setDerivedParameters()
      {
        if (Tg_ < TIME_CONSTANT_MINIMUM || TM_ < TIME_CONSTANT_MINIMUM)
        {
          static std::once_flag time_constant_warning_flag_;
          std::call_once(time_constant_warning_flag_,
                         &logTimeConstantWarning);
        }

        Tg_          = std::max(Tg_, TIME_CONSTANT_MINIMUM);
        TM_          = std::max(TM_, TIME_CONSTANT_MINIMUM);
        use_lvpl_    = ZERO<RealT>;
        bypass_lvpl_ = ONE<RealT>;
        if (sL_)
        {
          use_lvpl_    = ONE<RealT>;
          bypass_lvpl_ = ZERO<RealT>;
        }

        use_rqmax_ = ZERO<RealT>;
        use_rqmin_ = ZERO<RealT>;
        if (q0_ > ZERO<RealT> && Rqmax_ > ZERO<RealT>)
        {
          use_rqmax_ = ONE<RealT>;
        }
        else if (q0_ < ZERO<RealT> && Rqmin_ < ZERO<RealT>)
        {
          use_rqmin_ = ONE<RealT>;
        }

        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);
      }

      /**
       * @brief Convert a system-base per-unit value to the component base.
       *
       * @param[in] value Value on the system base.
       * @return Value on the REGCA component base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Regca<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * va_system_base_ / va_converter_base_;
      }

      /**
       * @brief Convert a component-base per-unit value to the system base.
       *
       * @param[in] value Value on the REGCA component base.
       * @return Value on the system base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Regca<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      /**
       * @brief Access the terminal-bus real voltage component.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Regca<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      /**
       * @brief Access the terminal-bus imaginary voltage component.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Regca<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

      /**
       * @brief Access the terminal-bus real-current accumulation target.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Regca<scalar_type, index_type>::Ir()
      {
        return bus_->Ir();
      }

      /**
       * @brief Access the terminal-bus imaginary-current accumulation target.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Regca<scalar_type, index_type>::Ii()
      {
        return bus_->Ii();
      }

      /**
       * @brief Calculate the initial HVRCM extra current.
       *
       * Solves the smooth HVRCM residual for \f$I_q^\mathrm{extra}\f$ given the
       * strictly positive voltage margin \f$V_\mathrm{hv}^{\max} - V_T\f$.
       *
       * @param[in] dv Strictly positive voltage margin.
       * @return Initial HVRCM extra current.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Regca<scalar_type, index_type>::initialHvrcmCurrent(
          scalar_type dv) const
      {
        static constexpr RealT log_two = std::numbers::ln2_v<RealT>;

        const ScalarT x = Math::MU<RealT> * dv;

        // Both branches evaluate log(1 - exp(-x)) and are algebraically
        // identical, so their values and derivatives agree at x = log(2).
        // The split only avoids cancellation for small x.
        if (x < log_two)
        {
          return -(log_two - HALF<RealT> * x
                   + std::log(std::sinh(HALF<RealT> * x)))
                 / Math::MU<RealT>;
        }

        return -std::log1p(-std::exp(-x)) / Math::MU<RealT>;
      }

      /**
       * @brief Set the LVPL release slope above the upper breakpoint.
       *
       * The finite slope is measured in p.u. current per p.u. voltage. The
       * unlimited characteristic is recovered as the slope approaches infinity.
       *
       * @param[in] KL Positive LVPL release slope.
       */
      template <typename scalar_type, typename index_type>
      void Regca<scalar_type, index_type>::setLvplGain(RealT KL)
      {
        KL_ = KL;
      }

      /**
       * @brief Read the required parameters from the model data.
       *
       * Missing keys and invalid value types are counted for verify(). Integer
       * JSON values are accepted for real parameters, and the LVPL switch accepts
       * a boolean or the integer values zero and one.
       *
       * @param[in] data REGCA parameter data.
       */
      template <typename scalar_type, typename index_type>
      void Regca<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        auto load_required_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Regca: missing required parameter '" << name << "'\n";
            ++parameter_error_count_;
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
            Log::error() << "Regca: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto load_required_switch = [&](auto key, bool& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Regca: missing required parameter '" << name << "'\n";
            ++parameter_error_count_;
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
          else
          {
            Log::error() << "Regca: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        load_required_real(Params::p0, p0_, "p0");
        load_required_real(Params::q0, q0_, "q0");
        load_required_real(Params::mva, mva_base_, "mva");
        load_required_real(Params::Tg, Tg_, "Tg");
        load_required_real(Params::TM, TM_, "TM");
        load_required_real(Params::Rqmax, Rqmax_, "Rqmax");
        load_required_real(Params::Rqmin, Rqmin_, "Rqmin");
        load_required_real(Params::Rpmax, Rpmax_, "Rpmax");
        load_required_switch(Params::sL, sL_, "sL");
        load_required_real(Params::IL1, IL1_, "IL1");
        load_required_real(Params::VL0, VL0_, "VL0");
        load_required_real(Params::VL1, VL1_, "VL1");
        load_required_real(Params::VA0, VA0_, "VA0");
        load_required_real(Params::VA1, VA1_, "VA1");
        load_required_real(Params::Vhvmax, Vhvmax_, "Vhvmax");

        setDerivedParameters();
      }

      /**
       * @brief Access the model monitor.
       *
       * @return Monitor for this model, or nullptr when constructed without data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Regca<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Bind monitorable variables to their internal states.
       *
       * Branch currents and powers are already on the system base and therefore
       * require no conversion before publication.
       */
      template <typename scalar_type, typename index_type>
      void Regca<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        auto index     = [](RegcaInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::ir, [this, index]
                      { return y_.getData()[index(RegcaInternalVariables::IR)]; });
        monitor_->set(Variable::ii, [this, index]
                      { return y_.getData()[index(RegcaInternalVariables::II)]; });
        monitor_->set(Variable::p, [this, index]
                      { return y_.getData()[index(RegcaInternalVariables::PBR)]; });
        monitor_->set(Variable::q, [this, index]
                      { return y_.getData()[index(RegcaInternalVariables::QBR)]; });
      }

      /**
       * @brief Set the component identifier assigned by the system model.
       *
       * @param[in] component_id Component identifier.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate model storage and connect assigned output signals.
       *
       * State, residual, tolerance, and interface buffers are sized; local index
       * maps are initialized; and assigned output signals are connected to the
       * states they publish. Repeated calls reuse the allocated model vectors.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.resize(2);
        wb_.setToZero();
        h_.resize(2);
        h_.setToZero();

        auto signal_size = static_cast<size_t>(RegcaExternalVariables::MAXIMUM);
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        auto* y = y_.getData();

        if (signals_.template isAssigned<RegcaInternalVariables::IR>())
        {
          signals_.template getSignalNode<RegcaInternalVariables::IR>()->set(
              &y[static_cast<size_t>(RegcaInternalVariables::IR)],
              &(this->getVariableIndex(static_cast<IdxT>(RegcaInternalVariables::IR))));
        }

        if (signals_.template isAssigned<RegcaInternalVariables::II>())
        {
          signals_.template getSignalNode<RegcaInternalVariables::II>()->set(
              &y[static_cast<size_t>(RegcaInternalVariables::II)],
              &(this->getVariableIndex(static_cast<IdxT>(RegcaInternalVariables::II))));
        }

        if (signals_.template isAssigned<RegcaInternalVariables::PBR>())
        {
          signals_.template getSignalNode<RegcaInternalVariables::PBR>()->set(
              &y[static_cast<size_t>(RegcaInternalVariables::PBR)],
              &(this->getVariableIndex(static_cast<IdxT>(RegcaInternalVariables::PBR))));
        }

        if (signals_.template isAssigned<RegcaInternalVariables::QBR>())
        {
          signals_.template getSignalNode<RegcaInternalVariables::QBR>()->set(
              &y[static_cast<size_t>(RegcaInternalVariables::QBR)],
              &(this->getVariableIndex(static_cast<IdxT>(RegcaInternalVariables::QBR))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the REGCA configuration.
       *
       * Applies the parameter conditions documented in the README, requires a
       * terminal bus, and checks that attached command ports have linked sources.
       * Operating-point admissibility is checked by initialize().
       *
       * @return Number of configuration errors, zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Regca: " << message << '\n';
            ret += 1;
          }
        };

        if (bus_ == nullptr)
        {
          Log::error() << "Regca: bus pointer is null\n";
          ret += 1;
        }

        check(mva_base_ > ZERO<RealT>, "mva must be positive");
        check(Rpmax_ >= ZERO<RealT>, "Rpmax must be non-negative");
        check(IL1_ >= ZERO<RealT>, "IL1 must be non-negative");
        check(KL_ > ZERO<RealT>, "LVPL release slope must be positive");
        check(ZERO<RealT> <= VL0_ && VL0_ < VL1_, "VL0/VL1 must satisfy 0 <= VL0 < VL1");
        check(ZERO<RealT> <= VA0_ && VA0_ < VA1_ && VA1_ < Vhvmax_,
              "VA0/VA1/Vhvmax must satisfy 0 <= VA0 < VA1 < Vhvmax");

        if (signals_.template isAttached<RegcaExternalVariables::IPCMD>())
        {
          if (!signals_.template isLinked<RegcaExternalVariables::IPCMD>())
          {
            Log::error() << "Regca: ipcmd signal attached with no linked source\n";
            ret += 1;
          }
        }

        if (signals_.template isAttached<RegcaExternalVariables::IQCMD>())
        {
          if (!signals_.template isLinked<RegcaExternalVariables::IQCMD>())
          {
            Log::error() << "Regca: iqcmd signal attached with no linked source\n";
            ret += 1;
          }
        }

        return ret;
      }

      /**
       * @brief Initialize REGCA from the power-flow operating point.
       *
       * Resolves the internal states and current-command setpoints from the
       * initialized terminal-bus voltage and system-base P0/Q0 injections, then
       * publishes the resolved commands to attached ports.
       *
       * @pre allocate() has completed, verify() reports no errors, and the
       *      terminal bus has been initialized.
       * @pre \f$V_{A1} \le V_{T,0} < V_\mathrm{hv}^{\max}\f$, and with LVPL
       *      enabled \f$I_{p,0} \le I_{L,0}\f$.
       * @post All internal derivatives are zero. Unattached command ports retain
       *       the resolved setpoints as constant commands.
       * @return Zero on success; nonzero when the operating point is rejected.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::initialize()
      {
        const auto VM      = static_cast<size_t>(RegcaInternalVariables::VM);
        const auto IQ      = static_cast<size_t>(RegcaInternalVariables::IQ);
        const auto IP      = static_cast<size_t>(RegcaInternalVariables::IP);
        const auto VT      = static_cast<size_t>(RegcaInternalVariables::VT);
        const auto IR      = static_cast<size_t>(RegcaInternalVariables::IR);
        const auto II      = static_cast<size_t>(RegcaInternalVariables::II);
        const auto IQEXTRA = static_cast<size_t>(RegcaInternalVariables::IQEXTRA);
        const auto IL      = static_cast<size_t>(RegcaInternalVariables::IL);
        const auto PBR     = static_cast<size_t>(RegcaInternalVariables::PBR);
        const auto QBR     = static_cast<size_t>(RegcaInternalVariables::QBR);
        auto*      y       = y_.getData();

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();
        const ScalarT vt = std::sqrt(vr * vr + vi * vi);

        if (vt <= ZERO<RealT>)
        {
          Log::error() << "Regca: terminal voltage magnitude must be positive at initialization\n";
          return 1;
        }
        if (vt < VA1_)
        {
          Log::error()
              << "Regca: terminal voltage magnitude must be at least VA1 at initialization\n";
          return 1;
        }
        if (vt >= Vhvmax_)
        {
          Log::error()
              << "Regca: terminal voltage magnitude must be below Vhvmax at initialization\n";
          return 1;
        }

        // P0 is a system-base power-flow injection. Resolve the component-base
        // active current through the LVACM network-interface gain.
        const ScalarT lvacm = Math::linseg(vt, VA0_, VA1_, ONE<RealT>);
        const ScalarT ip0   = toComponentBase(static_cast<ScalarT>(p0_) / vt) / lvacm;
        const ScalarT il0   = Math::linseg(vt, VL0_, VL1_, IL1_)
                            + KL_ * Math::ramp(vt - VL1_);

        if (sL_ && ip0 > il0)
        {
          Log::error() << "Regca: initial active current exceeds the LVPL ceiling\n";
          return 1;
        }
        const ScalarT ipcmd0 = ip0;

        // Solve the smooth HVRCM constraint and preserve the requested Q0. The
        // Vhvmax check above keeps the voltage margin strictly positive, so the
        // solve is always finite.
        const ScalarT dv       = Vhvmax_ - vt;
        const ScalarT iqextra0 = initialHvrcmCurrent(dv);
        const ScalarT qnet0    = toComponentBase(static_cast<ScalarT>(q0_) / vt);
        const ScalarT iqcmd0   = qnet0 + iqextra0;
        const ScalarT ir0      = (vi * qnet0 + vr * ip0 * lvacm) / vt;
        const ScalarT ii0      = (-vr * qnet0 + vi * ip0 * lvacm) / vt;

        y[VM]      = vt;
        y[VT]      = vt;
        y[IP]      = ip0;
        y[IQ]      = iqcmd0;
        y[IQEXTRA] = iqextra0;
        y[IL]      = il0;
        y[IR]      = toSystemBase(ir0);
        y[II]      = toSystemBase(ii0);
        y[PBR]     = vr * y[IR] + vi * y[II];
        y[QBR]     = vi * y[IR] - vr * y[II];

        ipcmd_set_ = toSystemBase(ipcmd0);
        iqcmd_set_ = toSystemBase(iqcmd0);

        // Publish the resolved system-base commands for downstream controller
        // initialization. Unattached ports retain these values as constant
        // commands.
        if (signals_.template isAttached<RegcaExternalVariables::IPCMD>())
        {
          signals_.template writeExternalVariable<RegcaExternalVariables::IPCMD>(ipcmd_set_);
        }
        if (signals_.template isAttached<RegcaExternalVariables::IQCMD>())
        {
          signals_.template writeExternalVariable<RegcaExternalVariables::IQCMD>(iqcmd_set_);
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief Tag VM, IQ, and IP as differential variables.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(RegcaInternalVariables::VM)] = true;
        tag_[static_cast<size_t>(RegcaInternalVariables::IQ)] = true;
        tag_[static_cast<size_t>(RegcaInternalVariables::IP)] = true;
        return 0;
      }

      /**
       * @brief Set every REGCA absolute tolerance to the solver relative tolerance.
       *
       * @param[in] rel_tol Solver relative tolerance.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Evaluate the REGCA internal residual equations.
       *
       * Evaluates the converter states and algebraic rows in
       * RegcaInternalVariables order. The body remains branchless so sparse
       * automatic differentiation sees a configuration-independent structure.
       *
       * @param[in] y Internal variables.
       * @param[in] yp Internal variable derivatives.
       * @param[in] wb Terminal-bus voltage components.
       * @param[in] ws Current-command signal values on the system base.
       * @param[out] f Internal residuals.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Regca<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto VM      = static_cast<size_t>(RegcaInternalVariables::VM);
        const auto IQ      = static_cast<size_t>(RegcaInternalVariables::IQ);
        const auto IP      = static_cast<size_t>(RegcaInternalVariables::IP);
        const auto VT      = static_cast<size_t>(RegcaInternalVariables::VT);
        const auto IR      = static_cast<size_t>(RegcaInternalVariables::IR);
        const auto II      = static_cast<size_t>(RegcaInternalVariables::II);
        const auto IQEXTRA = static_cast<size_t>(RegcaInternalVariables::IQEXTRA);
        const auto IL      = static_cast<size_t>(RegcaInternalVariables::IL);
        const auto PBR     = static_cast<size_t>(RegcaInternalVariables::PBR);
        const auto QBR     = static_cast<size_t>(RegcaInternalVariables::QBR);

        const auto IPCMD = static_cast<size_t>(RegcaExternalVariables::IPCMD);
        const auto IQCMD = static_cast<size_t>(RegcaExternalVariables::IQCMD);

        const ScalarT vm      = y[VM];
        const ScalarT iq      = y[IQ];
        const ScalarT ip      = y[IP];
        const ScalarT vt      = y[VT];
        const ScalarT ir      = y[IR];
        const ScalarT ii      = y[II];
        const ScalarT iqextra = y[IQEXTRA];
        const ScalarT il      = y[IL];
        const ScalarT pbr     = y[PBR];
        const ScalarT qbr     = y[QBR];

        const ScalarT vm_dot = yp[VM];
        const ScalarT iq_dot = yp[IQ];
        const ScalarT ip_dot = yp[IP];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT ipcmd = toComponentBase(ws[IPCMD]);
        const ScalarT iqcmd = toComponentBase(ws[IQCMD]);

        // Form the unconstrained current derivatives, then apply the REGCA
        // recovery rate limits in p.u./s.
        const ScalarT fq = (iqcmd - iq) / Tg_;
        const ScalarT fp = (ipcmd - ip) / Tg_;

        // At Q0 = 0 both corrections vanish, leaving fq unrestricted.
        const ScalarT iq_rate = fq + use_rqmax_ * (Math::min(fq, Rqmax_) - fq)
                                + use_rqmin_ * (Math::max(fq, Rqmin_) - fq);
        const ScalarT fp_limited = rrpwr(ip, fp, Rpmax_);

        // The LVPL ceiling IL = linseg(VM) moves with the sensed voltage; its
        // rate is the exact chain rule (inside() is the linseg slope mask
        // since ramp' = sigmoid), and a pinned Ip tracks the moving ceiling.
        const ScalarT vm_rate = (vt - vm) / TM_;
        const ScalarT il_rate = (IL1_ / (VL1_ - VL0_) * Math::inside(vm, VL0_, VL1_)
                                 + KL_ * Math::sigmoid(vm - VL1_))
                                * vm_rate;
        const ScalarT lvacm = Math::linseg(vt, VA0_, VA1_, ONE<RealT>);
        const ScalarT qnet  = iq - iqextra;

        f[VM] = -vm_dot + vm_rate;
        f[IQ] = -iq_dot + iq_rate;
        f[IP] = -ip_dot + bypass_lvpl_ * fp_limited
                + use_lvpl_ * awmax(ip, fp_limited, il, il_rate);
        f[VT]      = -vt * vt + vr * vr + vi * vi;
        f[IR]      = -toComponentBase(vt * ir) + vi * qnet + vr * ip * lvacm;
        f[II]      = -toComponentBase(vt * ii) - vr * qnet + vi * ip * lvacm;
        f[IQEXTRA] = -iqextra + Math::ramp(iqextra - (Vhvmax_ - vt));
        f[IL]      = -il + Math::linseg(vm, VL0_, VL1_, IL1_)
                + KL_ * Math::ramp(vm - VL1_);
        f[PBR] = -pbr + vr * ir + vi * ii;
        f[QBR] = -qbr + vi * ir - vr * ii;

        return 0;
      }

      /**
       * @brief Evaluate the current injected into the terminal bus.
       *
       * The branch-current states are already on the system base.
       *
       * @param[in] y Internal variables.
       * @param[in] yp Internal variable derivatives, unused.
       * @param[in] wb Terminal-bus voltage components, unused.
       * @param[out] h Current injected into the terminal bus.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int Regca<scalar_type, index_type>::evaluateBusResidual(
          const ScalarT*                  y,
          [[maybe_unused]] const ScalarT* yp,
          [[maybe_unused]] const ScalarT* wb,
          ScalarT*                        h)
      {
        const auto IR = static_cast<size_t>(RegcaInternalVariables::IR);
        const auto II = static_cast<size_t>(RegcaInternalVariables::II);

        h[0] = y[IR];
        h[1] = y[II];
        return 0;
      }

      /**
       * @brief Evaluate model residuals and accumulate the branch current.
       *
       * Refreshes the bus and signal interface buffers, evaluates the internal
       * and bus residuals, and accumulates the converter current into the terminal
       * bus. Unattached command ports use the setpoints latched by initialize().
       *
       * @pre The terminal-bus residual has been zeroed for this evaluation.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateResidual()
      {
        const auto IPCMD = static_cast<size_t>(RegcaExternalVariables::IPCMD);
        const auto IQCMD = static_cast<size_t>(RegcaExternalVariables::IQCMD);

        auto* ws = ws_.getData();

        ws[IPCMD] = ipcmd_set_;
        ws[IQCMD] = iqcmd_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<RegcaExternalVariables::IPCMD>())
        {
          ws[IPCMD] = signals_.template readExternalVariable<RegcaExternalVariables::IPCMD>();
          ws_indices_[IPCMD] =
              signals_.template readExternalVariableIndex<RegcaExternalVariables::IPCMD>();
        }

        if (signals_.template isAttached<RegcaExternalVariables::IQCMD>())
        {
          ws[IQCMD] = signals_.template readExternalVariable<RegcaExternalVariables::IQCMD>();
          ws_indices_[IQCMD] =
              signals_.template readExternalVariableIndex<RegcaExternalVariables::IQCMD>();
        }

        auto* wb = wb_.getData();
        wb[0]    = Vr();
        wb[1]    = Vi();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        auto*       h  = h_.getData();

        evaluateInternalResidual(y, yp, wb, ws, f);
        evaluateBusResidual(y, yp, wb, h);
        f_.setDataUpdated();

        Ir() += h[0];
        Ii() += h[1];
        bus_->getResidual().setDataUpdated();

        return 0;
      }
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

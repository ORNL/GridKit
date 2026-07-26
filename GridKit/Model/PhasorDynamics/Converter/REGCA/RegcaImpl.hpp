/**
 * @file RegcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REGCA phasor-dynamics converter model.
 */

#pragma once

#include <algorithm>
#include <cmath>
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

      template <typename scalar_type, typename index_type>
      Regca<scalar_type, index_type>::Regca(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
      }

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

      template <typename scalar_type, typename index_type>
      void Regca<scalar_type, index_type>::setDerivedParameters()
      {
        if (Tg_ < TIME_CONSTANT_MINIMUM || TM_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "Regca: Tg and TM below " << TIME_CONSTANT_MINIMUM
                         << " s are raised to that floor to keep the current-control "
                         << "and voltage-sensor lags well posed\n";
        }

        Tg_          = std::max(Tg_, TIME_CONSTANT_MINIMUM);
        TM_          = std::max(TM_, TIME_CONSTANT_MINIMUM);
        Mp_          = INACTIVE_RATE_LIMIT_FACTOR * Rpmax_;
        use_lvpl_    = ZERO<RealT>;
        bypass_lvpl_ = ONE<RealT>;
        if (sL_)
        {
          use_lvpl_    = ONE<RealT>;
          bypass_lvpl_ = ZERO<RealT>;
        }
        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);
      }

      template <typename scalar_type, typename index_type>
      scalar_type Regca<scalar_type, index_type>::lpTarget(scalar_type ip) const
      {
        return -Rpmax_ - (Mp_ - Rpmax_) * Math::sigmoid(ip);
      }

      template <typename scalar_type, typename index_type>
      scalar_type Regca<scalar_type, index_type>::upTarget(
          scalar_type ip,
          scalar_type il) const
      {
        const ScalarT sigma_ip = Math::sigmoid(ip);
        return Mp_ * (ONE<RealT> - sigma_ip)
               + Rpmax_ * sigma_ip * (bypass_lvpl_ + use_lvpl_ * Math::sigmoid(il - ip));
      }

      template <typename scalar_type, typename index_type>
      scalar_type Regca<scalar_type, index_type>::solveInitialHvrcmCurrent(
          scalar_type voltage_margin) const
      {
        static constexpr RealT log_two = static_cast<RealT>(0.6931471805599453);

        const ScalarT scaled_margin = Math::MU<RealT> * voltage_margin;

        // q = ramp(q - m) gives q = -log(1 - exp(-mu * m)) / mu.
        // Evaluate log(1 - exp(-x)) without cancellation.
        ScalarT log_one_minus_exp;
        if (scaled_margin < log_two)
        {
          log_one_minus_exp = log_two - HALF<RealT> * scaled_margin
                              + std::log(std::sinh(HALF<RealT> * scaled_margin));
        }
        else
        {
          log_one_minus_exp = std::log1p(-std::exp(-scaled_margin));
        }

        return -log_one_minus_exp / Math::MU<RealT>;
      }

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

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Regca<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

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

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

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

        wb_.assign(2, ScalarT{0});
        h_.assign(2, ScalarT{0});

        auto signal_size = static_cast<size_t>(RegcaExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
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
       * @brief Validate the REGCA configuration
       *
       * Checks required parameters, terminal-bus association, and attached
       * command ports. Operating-point-dependent constraints are part of
       * initialize()'s documented preconditions.
       *
       * @return int Number of configuration errors; zero when the configuration
       *             is valid.
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
        check(Rpmax_ > ZERO<RealT>, "Rpmax must be positive");
        check(Rqmin_ < ZERO<RealT> && ZERO<RealT> < Rqmax_, "Rqmin < 0 < Rqmax is required");
        check(IL1_ >= ZERO<RealT>, "IL1 must be non-negative");
        check(ZERO<RealT> <= VL0_ && VL0_ < VL1_, "VL0/VL1 must satisfy 0 <= VL0 < VL1");
        check(ZERO<RealT> <= VA0_ && VA0_ < VA1_, "VA0/VA1 must satisfy 0 <= VA0 < VA1");
        check(Vhvmax_ > ZERO<RealT>, "Vhvmax must be positive");

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
       * @brief Initialize REGCA from the power-flow operating point
       *
       * Resolves the internal states and current-command signals from the
       * initialized terminal-bus voltage and the system-base P0/Q0 injections.
       * All differential states are initialized with zero derivative.
       *
       * The initial terminal-voltage magnitude must be at least VA1 and below
       * Vhvmax.
       *
       * @pre allocate() has completed.
       * @pre verify() has reported no configuration errors.
       * @pre The terminal bus has been initialized.
       *
       * @return int 0 on success; nonzero when the terminal voltage is
       *             non-positive, below VA1, or not below Vhvmax.
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
        const auto LP      = static_cast<size_t>(RegcaInternalVariables::LP);
        const auto UP      = static_cast<size_t>(RegcaInternalVariables::UP);
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
        // active-current command through the LVACM network-interface gain.
        const ScalarT lvacm  = Math::linseg(vt, VA0_, VA1_, ONE<RealT>);
        const ScalarT ipcmd0 = toComponentBase(static_cast<ScalarT>(p0_) / vt) / lvacm;

        // Solve the smooth HVRCM constraint and preserve the requested Q0. The
        // Vhvmax check above keeps the voltage margin strictly positive, so the
        // solve is always finite.
        const ScalarT iqextra0 = solveInitialHvrcmCurrent(Vhvmax_ - vt);
        const ScalarT qnet0    = toComponentBase(static_cast<ScalarT>(q0_) / vt);
        const ScalarT iqcmd0   = qnet0 + iqextra0;
        const ScalarT ir0      = (vi * qnet0 + vr * ipcmd0 * lvacm) / vt;
        const ScalarT ii0      = (-vr * qnet0 + vi * ipcmd0 * lvacm) / vt;

        y[VM]      = vt;
        y[VT]      = vt;
        y[IP]      = ipcmd0;
        y[IQ]      = iqcmd0;
        y[IQEXTRA] = iqextra0;
        y[IL]      = Math::linseg(vt, VL0_, VL1_, IL1_);
        y[IR]      = toSystemBase(ir0);
        y[II]      = toSystemBase(ii0);
        y[LP]      = lpTarget(y[IP]);
        y[UP]      = upTarget(y[IP], y[IL]);
        y[PBR]     = vr * y[IR] + vi * y[II];
        y[QBR]     = vi * y[IR] - vr * y[II];

        // The WECC REGCA definition selects the reactive-current recovery
        // limiter from the sign of the initial reactive-power injection, not
        // from IQCMD after HVRCM compensation.
        ipcmd_set_    = toSystemBase(ipcmd0);
        iqcmd_set_    = toSystemBase(iqcmd0);
        iq_use_upper_ = ZERO<RealT>;
        iq_use_lower_ = ONE<RealT>;
        if (q0_ > ZERO<RealT>)
        {
          iq_use_upper_ = ONE<RealT>;
          iq_use_lower_ = ZERO<RealT>;
        }

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

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(RegcaInternalVariables::VM)] = true;
        tag_[static_cast<size_t>(RegcaInternalVariables::IQ)] = true;
        tag_[static_cast<size_t>(RegcaInternalVariables::IP)] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

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
        const auto LP      = static_cast<size_t>(RegcaInternalVariables::LP);
        const auto UP      = static_cast<size_t>(RegcaInternalVariables::UP);
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
        const ScalarT lp      = y[LP];
        const ScalarT up      = y[UP];
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

        const ScalarT iq_limited = iq_use_upper_ * Math::min(fq, Rqmax_)
                                   + iq_use_lower_ * Math::max(fq, Rqmin_);

        const ScalarT lvacm          = Math::linseg(vt, VA0_, VA1_, ONE<RealT>);
        const ScalarT qnet           = iq - iqextra;
        const ScalarT voltage_margin = Vhvmax_ - vt;

        f[VM]      = -vm_dot + (vt - vm) / TM_;
        f[IQ]      = -iq_dot + iq_limited;
        f[IP]      = -ip_dot + Math::clamp(fp, lp, up);
        f[VT]      = -vt * vt + vr * vr + vi * vi;
        f[IR]      = -toComponentBase(vt * ir) + vi * qnet + vr * ip * lvacm;
        f[II]      = -toComponentBase(vt * ii) - vr * qnet + vi * ip * lvacm;
        f[IQEXTRA] = -iqextra + Math::ramp(iqextra - voltage_margin);
        f[IL]      = -il + Math::linseg(vm, VL0_, VL1_, IL1_);
        f[LP]      = -lp + lpTarget(ip);
        f[UP]      = -up + upTarget(ip, il);
        f[PBR]     = -pbr + vr * ir + vi * ii;
        f[QBR]     = -qbr + vi * ir - vr * ii;

        return 0;
      }

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

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateResidual()
      {
        const auto IPCMD = static_cast<size_t>(RegcaExternalVariables::IPCMD);
        const auto IQCMD = static_cast<size_t>(RegcaExternalVariables::IQCMD);

        ws_[IPCMD] = ipcmd_set_;
        ws_[IQCMD] = iqcmd_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<RegcaExternalVariables::IPCMD>())
        {
          ws_[IPCMD] = signals_.template readExternalVariable<RegcaExternalVariables::IPCMD>();
          ws_indices_[IPCMD] =
              signals_.template readExternalVariableIndex<RegcaExternalVariables::IPCMD>();
        }

        if (signals_.template isAttached<RegcaExternalVariables::IQCMD>())
        {
          ws_[IQCMD] = signals_.template readExternalVariable<RegcaExternalVariables::IQCMD>();
          ws_indices_[IQCMD] =
              signals_.template readExternalVariableIndex<RegcaExternalVariables::IQCMD>();
        }

        wb_[0] = Vr();
        wb_[1] = Vi();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);
        evaluateBusResidual(y, yp, wb_.data(), h_.data());
        f_.setDataUpdated();

        Ir() += h_[0];
        Ii() += h_[1];
        bus_->getResidual().setDataUpdated();

        return 0;
      }
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

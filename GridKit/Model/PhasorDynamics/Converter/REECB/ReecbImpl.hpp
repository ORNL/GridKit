/**
 * @file ReecbImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REECB electrical-control model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/ReecbData.hpp>
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
      Reecb<scalar_type, index_type>::Reecb(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(ReecbInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::Reecb(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(ReecbInternalVariables::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::~Reecb()
      {
      }

      template <typename scalar_type, typename index_type>
      scalar_type& Reecb<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      template <typename scalar_type, typename index_type>
      scalar_type& Reecb<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::setDerivedParameters()
      {
        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);
        Trv_eff_           = std::max(Trv_, TIME_CONSTANT_MINIMUM);
        Tp_eff_            = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        pf_off_            = ONE<RealT> - PfFlag_;
        v_off_             = ONE<RealT> - VFlag_;
        q_off_             = ONE<RealT> - QFlag_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Reecb<scalar_type, index_type>::toComponentBase(
          scalar_type value) const
      {
        return value * va_system_base_ / va_converter_base_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Reecb<scalar_type, index_type>::toSystemBase(
          scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;
        Vref0_given_           = false;

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
            Log::error() << "Reecb: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto load_switch = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* bool_value = std::get_if<bool>(&value))
          {
            target = ZERO<RealT>;
            if (*bool_value)
            {
              target = ONE<RealT>;
            }
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value);
                   index_value && (*index_value == 0 || *index_value == 1))
          {
            target = static_cast<RealT>(*index_value);
          }
          else if (const auto* real_value = std::get_if<RealT>(&value);
                   real_value && (*real_value == ZERO<RealT> || *real_value == ONE<RealT>) )
          {
            target = *real_value;
          }
          else
          {
            Log::error() << "Reecb: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        if (!data.parameters.contains(Params::mva))
        {
          Log::error() << "Reecb: missing required parameter 'mva'\n";
          ++parameter_error_count_;
        }
        load_real(Params::mva, mva_base_, "mva");
        load_switch(Params::PfFlag, PfFlag_, "PfFlag");
        load_switch(Params::VFlag, VFlag_, "VFlag");
        load_switch(Params::QFlag, QFlag_, "QFlag");
        load_switch(Params::Pqflag, Pqflag_, "Pqflag");
        load_real(Params::Trv, Trv_, "Trv");
        load_real(Params::Tp, Tp_, "Tp");
        if (data.parameters.contains(Params::Vref0))
        {
          load_real(Params::Vref0, Vref0_, "Vref0");
          Vref0_given_ = true;
        }
        load_real(Params::Vdip, Vdip_, "Vdip");
        load_real(Params::Vup, Vup_, "Vup");
        load_real(Params::dbd1, dbd1_, "dbd1");
        load_real(Params::dbd2, dbd2_, "dbd2");
        load_real(Params::kqv, kqv_, "kqv");
        load_real(Params::Iql1, Iql1_, "Iql1");
        load_real(Params::Iqh1, Iqh1_, "Iqh1");
        load_real(Params::Qmax, Qmax_, "Qmax");
        load_real(Params::Qmin, Qmin_, "Qmin");
        load_real(Params::Kqp, Kqp_, "Kqp");
        load_real(Params::Kqi, Kqi_, "Kqi");
        load_real(Params::Vmax, Vmax_, "Vmax");
        load_real(Params::Vmin, Vmin_, "Vmin");
        load_real(Params::Kvp, Kvp_, "Kvp");
        load_real(Params::Kvi, Kvi_, "Kvi");
        load_real(Params::Tiq, Tiq_, "Tiq");
        load_real(Params::Tpord, Tpord_, "Tpord");
        load_real(Params::dPmax, dPmax_, "dPmax");
        load_real(Params::dPmin, dPmin_, "dPmin");
        load_real(Params::Pmax, Pmax_, "Pmax");
        load_real(Params::Pmin, Pmin_, "Pmin");
        load_real(Params::Imax, Imax_, "Imax");
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Reecb<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        auto index     = [](ReecbInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::iqcmd, [this, index]
                      { return y_.getData()[index(ReecbInternalVariables::IQCMD)]; });
        monitor_->set(Variable::ipcmd, [this, index]
                      { return y_.getData()[index(ReecbInternalVariables::IPCMD)]; });
        monitor_->set(Variable::vmeas, [this, index]
                      { return y_.getData()[index(ReecbInternalVariables::VMEAS)]; });
        monitor_->set(Variable::pmeas, [this, index]
                      { return y_.getData()[index(ReecbInternalVariables::PMEAS)]; });
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::allocate()
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

        auto signal_size = static_cast<size_t>(ReecbExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<ReecbInternalVariables::IQCMD>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<ReecbInternalVariables::IQCMD>()->set(
              &y[static_cast<size_t>(ReecbInternalVariables::IQCMD)],
              &(this->getVariableIndex(static_cast<IdxT>(ReecbInternalVariables::IQCMD))));
        }

        if (signals_.template isAssigned<ReecbInternalVariables::IPCMD>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<ReecbInternalVariables::IPCMD>()->set(
              &y[static_cast<size_t>(ReecbInternalVariables::IPCMD)],
              &(this->getVariableIndex(static_cast<IdxT>(ReecbInternalVariables::IPCMD))));
        }

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Reecb: " << message << '\n';
            ret += 1;
          }
        };

        if (bus_ == nullptr)
        {
          Log::error() << "Reecb: bus pointer is null\n";
          ret += 1;
        }

        check(mva_base_ > ZERO<RealT>, "mva must be positive");
        check(va_converter_base_ > ZERO<RealT>, "converter VA base must be positive");
        check(PfFlag_ == ZERO<RealT> || PfFlag_ == ONE<RealT>, "PfFlag must be 0 or 1");
        check(VFlag_ == ZERO<RealT> || VFlag_ == ONE<RealT>, "VFlag must be 0 or 1");
        check(QFlag_ == ZERO<RealT> || QFlag_ == ONE<RealT>, "QFlag must be 0 or 1");
        check(Pqflag_ == ZERO<RealT> || Pqflag_ == ONE<RealT>, "Pqflag must be 0 or 1");
        check(Trv_ >= ZERO<RealT>, "Trv must be non-negative");
        check(Tp_ >= ZERO<RealT>, "Tp must be non-negative");
        check(Vdip_ < Vup_, "Vdip must be less than Vup");
        check(dbd1_ <= ZERO<RealT> && ZERO<RealT> <= dbd2_, "dbd1 <= 0 <= dbd2 is required");
        check(Iql1_ <= Iqh1_, "Iql1 must be less than or equal to Iqh1");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax");
        check(Tiq_ > ZERO<RealT>, "Tiq must be positive");
        check(Tpord_ > ZERO<RealT>, "Tpord must be positive");
        check(dPmin_ < ZERO<RealT> && ZERO<RealT> < dPmax_, "dPmin < 0 < dPmax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");
        check(Imax_ >= ZERO<RealT>, "Imax must be non-negative");

        auto check_attached_signal = [&](bool attached, bool linked, const char* name)
        {
          if (attached && !linked)
          {
            Log::error() << "Reecb: " << name << " signal attached with no linked variable\n";
            ret += 1;
          }
        };

        check_attached_signal(
            signals_.template isAttached<ReecbExternalVariables::PE>(),
            signals_.template isAttached<ReecbExternalVariables::PE>()
                && signals_.template isLinked<ReecbExternalVariables::PE>(),
            "pe");
        check_attached_signal(
            signals_.template isAttached<ReecbExternalVariables::QGEN>(),
            signals_.template isAttached<ReecbExternalVariables::QGEN>()
                && signals_.template isLinked<ReecbExternalVariables::QGEN>(),
            "qgen");
        check_attached_signal(
            signals_.template isAttached<ReecbExternalVariables::QEXT>(),
            signals_.template isAttached<ReecbExternalVariables::QEXT>()
                && signals_.template isLinked<ReecbExternalVariables::QEXT>(),
            "qext");
        check_attached_signal(
            signals_.template isAttached<ReecbExternalVariables::PFAREF>(),
            signals_.template isAttached<ReecbExternalVariables::PFAREF>()
                && signals_.template isLinked<ReecbExternalVariables::PFAREF>(),
            "pfaref");
        check_attached_signal(
            signals_.template isAttached<ReecbExternalVariables::PREF>(),
            signals_.template isAttached<ReecbExternalVariables::PREF>()
                && signals_.template isLinked<ReecbExternalVariables::PREF>(),
            "pref");

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::initialize()
      {
        if (parameter_error_count_ > 0 || verify() > 0)
        {
          Log::error() << "Reecb: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        const auto VMEAS     = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS     = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ      = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV      = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV        = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD      = static_cast<size_t>(ReecbInternalVariables::PORD);
        const auto VT        = static_cast<size_t>(ReecbInternalVariables::VT);
        const auto VMEASSAFE = static_cast<size_t>(ReecbInternalVariables::VMEASSAFE);
        const auto SDIP      = static_cast<size_t>(ReecbInternalVariables::SDIP);
        const auto VERR      = static_cast<size_t>(ReecbInternalVariables::VERR);
        const auto IQV       = static_cast<size_t>(ReecbInternalVariables::IQV);
        const auto QREF      = static_cast<size_t>(ReecbInternalVariables::QREF);
        const auto EQ        = static_cast<size_t>(ReecbInternalVariables::EQ);
        const auto VPIQ      = static_cast<size_t>(ReecbInternalVariables::VPIQ);
        const auto EPIV      = static_cast<size_t>(ReecbInternalVariables::EPIV);
        const auto FPORD     = static_cast<size_t>(ReecbInternalVariables::FPORD);
        const auto RPORD     = static_cast<size_t>(ReecbInternalVariables::RPORD);
        const auto IQCIRC    = static_cast<size_t>(ReecbInternalVariables::IQCIRC);
        const auto IPCIRC    = static_cast<size_t>(ReecbInternalVariables::IPCIRC);
        const auto IQMAX     = static_cast<size_t>(ReecbInternalVariables::IQMAX);
        const auto IPMAX     = static_cast<size_t>(ReecbInternalVariables::IPMAX);
        const auto IQBASE    = static_cast<size_t>(ReecbInternalVariables::IQBASE);
        const auto IQRAW     = static_cast<size_t>(ReecbInternalVariables::IQRAW);
        const auto IQCMD     = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD     = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();

        y[VT] = std::sqrt(vr * vr + vi * vi);
        if (!Vref0_given_)
        {
          Vref0_ = static_cast<RealT>(y[VT]);
        }

        y[VMEAS]     = y[VT];
        y[VMEASSAFE] = Math::max(y[VMEAS], VMEAS_MINIMUM);

        const ScalarT ipcmd0 = toComponentBase(y[IPCMD]);
        const ScalarT iqcmd0 = toComponentBase(y[IQCMD]);
        ScalarT       pe0    = ipcmd0 * y[VMEASSAFE];
        ScalarT       qgen0  = iqcmd0 * y[VMEASSAFE];

        const ScalarT qext0 = qgen0;
        const ScalarT pref0 = Math::clamp(pe0, Pmin_, Pmax_);

        pe_set_     = toSystemBase(pe0);
        qgen_set_   = toSystemBase(qgen0);
        qext_set_   = toSystemBase(qext0);
        pfaref_set_ = std::abs(static_cast<RealT>(pe0)) > INIT_TOL ? static_cast<ScalarT>(std::atan(static_cast<RealT>(qgen0 / pe0))) : static_cast<ScalarT>(ZERO<RealT>);
        pref_set_   = toSystemBase(pref0);

        if (signals_.template isAttached<ReecbExternalVariables::PE>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::PE>(pe_set_);
        }
        if (signals_.template isAttached<ReecbExternalVariables::QGEN>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::QGEN>(qgen_set_);
        }
        if (signals_.template isAttached<ReecbExternalVariables::QEXT>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::QEXT>(qext_set_);
        }
        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::PFAREF>(pfaref_set_);
        }
        if (signals_.template isAttached<ReecbExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::PREF>(pref_set_);
        }

        y[PMEAS] = pe0;
        y[SDIP]  = Math::inside(y[VT], Vdip_, Vup_);
        y[VERR]  = Math::deadband2(Vref0_ - y[VMEAS], dbd1_, dbd2_);
        y[IQV]   = Math::clamp(kqv_ * y[VERR], Iql1_, Iqh1_);
        y[QREF]  = PfFlag_ * y[PMEAS] * std::tan(pfaref_set_) + pf_off_ * qext0;
        y[EQ]    = Math::clamp(y[QREF], Qmin_, Qmax_) - qgen0;
        y[QV]    = y[QREF] / y[VMEASSAFE];
        y[PORD]  = pref0;
        y[FPORD] = ZERO<RealT>;
        y[RPORD] = ZERO<RealT>;

        auto awinit = [](const ScalarT target, const ScalarT rate, const ScalarT lower, const ScalarT upper) -> ScalarT
        {
          if (std::abs(static_cast<RealT>(rate)) <= INIT_TOL)
          {
            return target;
          }
          return rate > ZERO<RealT> ? upper + static_cast<ScalarT>(SAT_MARGIN) : lower - static_cast<ScalarT>(SAT_MARGIN);
        };

        const ScalarT vpiq_arg = awinit(VFlag_ * y[VMEAS] + v_off_ * y[QREF], Kqi_ * y[EQ], static_cast<ScalarT>(Vmin_), static_cast<ScalarT>(Vmax_));
        y[VPIQ]                = Math::clamp(vpiq_arg, Vmin_, Vmax_);
        y[EPIV]                = VFlag_ * y[VPIQ] + v_off_ * y[QREF] - y[VMEAS];
        y[XPIQ]                = vpiq_arg - Kqp_ * y[EQ];

        const ScalarT iqbase_target = qgen0 / y[VMEASSAFE];
        const ScalarT ip_star       = y[PORD] / y[VMEASSAFE];

        auto initializeReactiveBase = [&]()
        {
          const ScalarT piv_arg = awinit(iqbase_target, Kvi_ * y[EPIV], -y[IQMAX], y[IQMAX]);
          y[IQBASE]             = Math::clamp(piv_arg, -y[IQMAX], y[IQMAX]);
          y[XPIV]               = piv_arg - Kvp_ * y[EPIV];
        };

        auto initializeReactiveCommand = [&]()
        {
          y[IQRAW] = QFlag_ * y[IQBASE] + q_off_ * y[QV] + (ONE<RealT> - y[SDIP]) * y[IQV];
          y[IQCMD] = toSystemBase(Math::clamp(y[IQRAW], -y[IQMAX], y[IQMAX]));
        };

        if (Pqflag_ == ZERO<RealT>)
        {
          y[IQCIRC] = Imax_;
          y[IQMAX]  = Imax_;
          initializeReactiveBase();
          initializeReactiveCommand();

          const ScalarT iqcmd       = toComponentBase(y[IQCMD]);
          const ScalarT ip_radicand = Imax_ * Imax_ - iqcmd * iqcmd;
          if (static_cast<RealT>(ip_radicand) < ZERO<RealT>)
          {
            Log::error() << "Reecb: initial active-current circle radicand is negative\n";
            return 1;
          }
          y[IPCIRC] = std::sqrt(ip_radicand);
          y[IPMAX]  = y[IPCIRC];
          y[IPCMD]  = toSystemBase(Math::clamp(ip_star, ZERO<RealT>, y[IPMAX]));
        }
        else
        {
          y[IPCIRC] = Imax_;
          y[IPMAX]  = Imax_;
          y[IPCMD]  = toSystemBase(Math::clamp(ip_star, ZERO<RealT>, y[IPMAX]));

          const ScalarT ipcmd       = toComponentBase(y[IPCMD]);
          const ScalarT iq_radicand = Imax_ * Imax_ - ipcmd * ipcmd;
          if (static_cast<RealT>(iq_radicand) < ZERO<RealT>)
          {
            Log::error() << "Reecb: initial reactive-current circle radicand is negative\n";
            return 1;
          }
          y[IQCIRC] = std::sqrt(iq_radicand);
          y[IQMAX]  = y[IQCIRC];
          initializeReactiveBase();
          initializeReactiveCommand();
        }

        for (IdxT i = 0; i < yp_.getSize(); ++i)
        {
          yp[i] = ZERO<RealT>;
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(ReecbInternalVariables::VMEAS)] = true;
        tag_[static_cast<size_t>(ReecbInternalVariables::PMEAS)] = true;
        tag_[static_cast<size_t>(ReecbInternalVariables::XPIQ)]  = true;
        tag_[static_cast<size_t>(ReecbInternalVariables::XPIV)]  = true;
        tag_[static_cast<size_t>(ReecbInternalVariables::QV)]    = true;
        tag_[static_cast<size_t>(ReecbInternalVariables::PORD)]  = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Reecb<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto VMEAS     = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS     = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ      = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV      = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV        = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD      = static_cast<size_t>(ReecbInternalVariables::PORD);
        const auto VT        = static_cast<size_t>(ReecbInternalVariables::VT);
        const auto VMEASSAFE = static_cast<size_t>(ReecbInternalVariables::VMEASSAFE);
        const auto SDIP      = static_cast<size_t>(ReecbInternalVariables::SDIP);
        const auto VERR      = static_cast<size_t>(ReecbInternalVariables::VERR);
        const auto IQV       = static_cast<size_t>(ReecbInternalVariables::IQV);
        const auto QREF      = static_cast<size_t>(ReecbInternalVariables::QREF);
        const auto EQ        = static_cast<size_t>(ReecbInternalVariables::EQ);
        const auto VPIQ      = static_cast<size_t>(ReecbInternalVariables::VPIQ);
        const auto EPIV      = static_cast<size_t>(ReecbInternalVariables::EPIV);
        const auto FPORD     = static_cast<size_t>(ReecbInternalVariables::FPORD);
        const auto RPORD     = static_cast<size_t>(ReecbInternalVariables::RPORD);
        const auto IQCIRC    = static_cast<size_t>(ReecbInternalVariables::IQCIRC);
        const auto IPCIRC    = static_cast<size_t>(ReecbInternalVariables::IPCIRC);
        const auto IQMAX     = static_cast<size_t>(ReecbInternalVariables::IQMAX);
        const auto IPMAX     = static_cast<size_t>(ReecbInternalVariables::IPMAX);
        const auto IQBASE    = static_cast<size_t>(ReecbInternalVariables::IQBASE);
        const auto IQRAW     = static_cast<size_t>(ReecbInternalVariables::IQRAW);
        const auto IQCMD     = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD     = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        const auto PE     = static_cast<size_t>(ReecbExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecbExternalVariables::QGEN);
        const auto QEXT   = static_cast<size_t>(ReecbExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecbExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecbExternalVariables::PREF);

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT pe     = toComponentBase(ws[PE]);
        const ScalarT qgen   = toComponentBase(ws[QGEN]);
        const ScalarT qext   = toComponentBase(ws[QEXT]);
        const ScalarT pfaref = ws[PFAREF];
        const ScalarT pref   = toComponentBase(ws[PREF]);
        const ScalarT iqcmd  = toComponentBase(y[IQCMD]);
        const ScalarT ipcmd  = toComponentBase(y[IPCMD]);

        f[VMEAS]     = -yp[VMEAS] + (y[VT] - y[VMEAS]) / Trv_eff_;
        f[PMEAS]     = -yp[PMEAS] + (pe - y[PMEAS]) / Tp_eff_;
        f[XPIQ]      = -yp[XPIQ] + y[SDIP] * Math::antiwindup(Kqp_ * y[EQ] + y[XPIQ], Kqi_ * y[EQ], Vmin_, Vmax_);
        f[XPIV]      = -yp[XPIV] + y[SDIP] * Math::antiwindup(Kvp_ * y[EPIV] + y[XPIV], Kvi_ * y[EPIV], -y[IQMAX], y[IQMAX]);
        f[QV]        = -yp[QV] + y[SDIP] * (y[QREF] / y[VMEASSAFE] - y[QV]) / Tiq_;
        f[PORD]      = -yp[PORD] + y[SDIP] * Math::antiwindup(y[PORD], y[RPORD], Pmin_, Pmax_);
        f[VT]        = -y[VT] * y[VT] + vr * vr + vi * vi;
        f[VMEASSAFE] = -y[VMEASSAFE] + Math::max(y[VMEAS], VMEAS_MINIMUM);
        f[SDIP]      = -y[SDIP] + Math::inside(y[VT], Vdip_, Vup_);
        f[VERR]      = -y[VERR] + Math::deadband2(Vref0_ - y[VMEAS], dbd1_, dbd2_);
        f[IQV]       = -y[IQV] + Math::clamp(kqv_ * y[VERR], Iql1_, Iqh1_);
        f[QREF]      = -y[QREF] + PfFlag_ * y[PMEAS] * std::tan(pfaref) + pf_off_ * qext;
        f[EQ]        = -y[EQ] + Math::clamp(y[QREF], Qmin_, Qmax_) - qgen;
        f[VPIQ]      = -y[VPIQ] + Math::clamp(Kqp_ * y[EQ] + y[XPIQ], Vmin_, Vmax_);
        f[EPIV]      = -y[EPIV] + VFlag_ * y[VPIQ] + v_off_ * y[QREF] - y[VMEAS];
        f[FPORD]     = -y[FPORD] + (pref - y[PORD]) / Tpord_;
        f[RPORD]     = -y[RPORD] + Math::clamp(y[FPORD], dPmin_, dPmax_);
        f[IQCIRC]    = -y[IQCIRC] * y[IQCIRC] + Imax_ * Imax_ - Pqflag_ * ipcmd * ipcmd;
        f[IPCIRC]    = -y[IPCIRC] * y[IPCIRC] + Imax_ * Imax_ - (ONE<RealT> - Pqflag_) * iqcmd * iqcmd;
        f[IQMAX]     = -y[IQMAX] + (ONE<RealT> - Pqflag_) * Imax_ + Pqflag_ * y[IQCIRC];
        f[IPMAX]     = -y[IPMAX] + Pqflag_ * Imax_ + (ONE<RealT> - Pqflag_) * y[IPCIRC];
        f[IQBASE]    = -y[IQBASE] + Math::clamp(Kvp_ * y[EPIV] + y[XPIV], -y[IQMAX], y[IQMAX]);
        f[IQRAW]     = -y[IQRAW] + QFlag_ * y[IQBASE] + q_off_ * y[QV] + (ONE<RealT> - y[SDIP]) * y[IQV];
        f[IQCMD]     = -y[IQCMD] + toSystemBase(Math::clamp(y[IQRAW], -y[IQMAX], y[IQMAX]));
        f[IPCMD]     = -y[IPCMD] + toSystemBase(Math::clamp(y[PORD] / y[VMEASSAFE], ZERO<RealT>, y[IPMAX]));

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateResidual()
      {
        const auto PE     = static_cast<size_t>(ReecbExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecbExternalVariables::QGEN);
        const auto QEXT   = static_cast<size_t>(ReecbExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecbExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecbExternalVariables::PREF);

        ws_[PE]     = pe_set_;
        ws_[QGEN]   = qgen_set_;
        ws_[QEXT]   = qext_set_;
        ws_[PFAREF] = pfaref_set_;
        ws_[PREF]   = pref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<ReecbExternalVariables::PE>())
        {
          ws_[PE]         = signals_.template readExternalVariable<ReecbExternalVariables::PE>();
          ws_indices_[PE] = signals_.template readExternalVariableIndex<ReecbExternalVariables::PE>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QGEN>())
        {
          ws_[QGEN]         = signals_.template readExternalVariable<ReecbExternalVariables::QGEN>();
          ws_indices_[QGEN] = signals_.template readExternalVariableIndex<ReecbExternalVariables::QGEN>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QEXT>())
        {
          ws_[QEXT]         = signals_.template readExternalVariable<ReecbExternalVariables::QEXT>();
          ws_indices_[QEXT] = signals_.template readExternalVariableIndex<ReecbExternalVariables::QEXT>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>())
        {
          ws_[PFAREF]         = signals_.template readExternalVariable<ReecbExternalVariables::PFAREF>();
          ws_indices_[PFAREF] = signals_.template readExternalVariableIndex<ReecbExternalVariables::PFAREF>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PREF>())
        {
          ws_[PREF]         = signals_.template readExternalVariable<ReecbExternalVariables::PREF>();
          ws_indices_[PREF] = signals_.template readExternalVariableIndex<ReecbExternalVariables::PREF>();
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
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

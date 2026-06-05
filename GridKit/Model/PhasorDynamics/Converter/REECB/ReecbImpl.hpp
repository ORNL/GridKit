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

      template <class ScalarT, typename IdxT>
      Reecb<ScalarT, IdxT>::Reecb(bus_type* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(ReecbInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Reecb<ScalarT, IdxT>::Reecb(bus_type* bus, const model_data_type& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(ReecbInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Reecb<ScalarT, IdxT>::~Reecb()
      {
      }

      template <class ScalarT, typename IdxT>
      ScalarT& Reecb<ScalarT, IdxT>::Vr()
      {
        return bus_->Vr();
      }

      template <class ScalarT, typename IdxT>
      ScalarT& Reecb<ScalarT, IdxT>::Vi()
      {
        return bus_->Vi();
      }

      template <class ScalarT, typename IdxT>
      void Reecb<ScalarT, IdxT>::setDerivedParameters()
      {
        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);
        pf_off_            = ONE<RealT> - PfFlag_;
        v_off_             = ONE<RealT> - VFlag_;
        q_off_             = ONE<RealT> - QFlag_;
      }

      template <class ScalarT, typename IdxT>
      ScalarT Reecb<ScalarT, IdxT>::toComponentBase(ScalarT value) const
      {
        return value * va_system_base_ / va_converter_base_;
      }

      template <class ScalarT, typename IdxT>
      void Reecb<ScalarT, IdxT>::initModelParams(const model_data_type& data)
      {
        using Params = typename model_data_type::Parameters;

        parameter_error_count_ = 0;
        Vref0_given_           = false;

        auto loadRequiredReal = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Reecb: missing required parameter '" << name << "'\n";
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
            Log::error() << "Reecb: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto loadRequiredSwitch = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Reecb: missing required parameter '" << name << "'\n";
            ++parameter_error_count_;
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

        loadRequiredReal(Params::mva, mva_base_, "mva");
        loadRequiredSwitch(Params::PfFlag, PfFlag_, "PfFlag");
        loadRequiredSwitch(Params::VFlag, VFlag_, "VFlag");
        loadRequiredSwitch(Params::QFlag, QFlag_, "QFlag");
        loadRequiredSwitch(Params::Pqflag, Pqflag_, "Pqflag");
        loadRequiredReal(Params::Trv, Trv_, "Trv");
        loadRequiredReal(Params::Tp, Tp_, "Tp");
        if (data.parameters.contains(Params::Vref0))
        {
          loadRequiredReal(Params::Vref0, Vref0_, "Vref0");
          Vref0_given_ = true;
        }
        loadRequiredReal(Params::Vdip, Vdip_, "Vdip");
        loadRequiredReal(Params::Vup, Vup_, "Vup");
        loadRequiredReal(Params::dbd1, dbd1_, "dbd1");
        loadRequiredReal(Params::dbd2, dbd2_, "dbd2");
        loadRequiredReal(Params::kqv, kqv_, "kqv");
        loadRequiredReal(Params::Iql1, Iql1_, "Iql1");
        loadRequiredReal(Params::Iqh1, Iqh1_, "Iqh1");
        loadRequiredReal(Params::Qmax, Qmax_, "Qmax");
        loadRequiredReal(Params::Qmin, Qmin_, "Qmin");
        loadRequiredReal(Params::Kqp, Kqp_, "Kqp");
        loadRequiredReal(Params::Kqi, Kqi_, "Kqi");
        loadRequiredReal(Params::Vmax, Vmax_, "Vmax");
        loadRequiredReal(Params::Vmin, Vmin_, "Vmin");
        loadRequiredReal(Params::Kvp, Kvp_, "Kvp");
        loadRequiredReal(Params::Kvi, Kvi_, "Kvi");
        loadRequiredReal(Params::Tiq, Tiq_, "Tiq");
        loadRequiredReal(Params::Tpord, Tpord_, "Tpord");
        loadRequiredReal(Params::dPmax, dPmax_, "dPmax");
        loadRequiredReal(Params::dPmin, dPmin_, "dPmin");
        loadRequiredReal(Params::Pmax, Pmax_, "Pmax");
        loadRequiredReal(Params::Pmin, Pmin_, "Pmin");
        loadRequiredReal(Params::Imax, Imax_, "Imax");
        setDerivedParameters();
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* Reecb<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void Reecb<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;
        auto index     = [](ReecbInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::iqcmd, [this, index]
                      { return y_[index(ReecbInternalVariables::IQCMD)]; });
        monitor_->set(Variable::ipcmd, [this, index]
                      { return y_[index(ReecbInternalVariables::IPCMD)]; });
        monitor_->set(Variable::vmeas, [this, index]
                      { return y_[index(ReecbInternalVariables::VMEAS)]; });
        monitor_->set(Variable::pmeas, [this, index]
                      { return y_[index(ReecbInternalVariables::PMEAS)]; });
        monitor_->set(Variable::piq, [this, index]
                      { return y_[index(ReecbInternalVariables::XPIQ)]; });
        monitor_->set(Variable::piv, [this, index]
                      { return y_[index(ReecbInternalVariables::XPIV)]; });
        monitor_->set(Variable::qv, [this, index]
                      { return y_[index(ReecbInternalVariables::QV)]; });
        monitor_->set(Variable::pord, [this, index]
                      { return y_[index(ReecbInternalVariables::PORD)]; });
        monitor_->set(Variable::qref, [this, index]
                      { return y_[index(ReecbInternalVariables::QREF)]; });
        monitor_->set(Variable::sdip, [this, index]
                      { return y_[index(ReecbInternalVariables::SDIP)]; });
        monitor_->set(Variable::iqmax, [this, index]
                      { return y_[index(ReecbInternalVariables::IQMAX)]; });
        monitor_->set(Variable::ipmax, [this, index]
                      { return y_[index(ReecbInternalVariables::IPMAX)]; });
        monitor_->set(Variable::iqv, [this, index]
                      { return y_[index(ReecbInternalVariables::IQV)]; });
        monitor_->set(Variable::vqctrl, [this, index]
                      { return y_[index(ReecbInternalVariables::VPIQ)]; });
        monitor_->set(Variable::iqbase, [this, index]
                      { return y_[index(ReecbInternalVariables::IQBASE)]; });
      }

      template <class ScalarT, typename IdxT>
      int Reecb<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Reecb<ScalarT, IdxT>::allocate()
      {
        size_     = static_cast<IdxT>(ReecbInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        f_.assign(size, ScalarT{0});
        y_.assign(size, ScalarT{0});
        yp_.assign(size, ScalarT{0});
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
          signals_.template getSignalNode<ReecbInternalVariables::IQCMD>()->set(
              &y_[static_cast<size_t>(ReecbInternalVariables::IQCMD)],
              &(this->getVariableIndex(static_cast<IdxT>(ReecbInternalVariables::IQCMD))));
        }

        if (signals_.template isAssigned<ReecbInternalVariables::IPCMD>())
        {
          signals_.template getSignalNode<ReecbInternalVariables::IPCMD>()->set(
              &y_[static_cast<size_t>(ReecbInternalVariables::IPCMD)],
              &(this->getVariableIndex(static_cast<IdxT>(ReecbInternalVariables::IPCMD))));
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Reecb<ScalarT, IdxT>::verify() const
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

        const bool has_pe   = signals_.template isAttached<ReecbExternalVariables::PE>();
        const bool has_qgen = signals_.template isAttached<ReecbExternalVariables::QGEN>();

        if (has_pe != has_qgen)
        {
          Log::error() << "Reecb: pe and qgen must be connected together\n";
          ret += 1;
        }

        if (signals_.template isAttached<ReecbExternalVariables::PE>()
            && !signals_.template isLinked<ReecbExternalVariables::PE>())
        {
          Log::error() << "Reecb: pe signal attached with no linked source\n";
          ret += 1;
        }

        if (signals_.template isAttached<ReecbExternalVariables::QGEN>()
            && !signals_.template isLinked<ReecbExternalVariables::QGEN>())
        {
          Log::error() << "Reecb: qgen signal attached with no linked source\n";
          ret += 1;
        }

        if (signals_.template isAttached<ReecbExternalVariables::QEXT>()
            && !signals_.template isLinked<ReecbExternalVariables::QEXT>())
        {
          Log::error() << "Reecb: qext signal attached with no linked source\n";
          ret += 1;
        }

        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>()
            && !signals_.template isLinked<ReecbExternalVariables::PFAREF>())
        {
          Log::error() << "Reecb: pfaref signal attached with no linked source\n";
          ret += 1;
        }

        if (signals_.template isAttached<ReecbExternalVariables::PREF>()
            && !signals_.template isLinked<ReecbExternalVariables::PREF>())
        {
          Log::error() << "Reecb: pref signal attached with no linked source\n";
          ret += 1;
        }

        return ret;
      }

      template <class ScalarT, typename IdxT>
      int Reecb<ScalarT, IdxT>::initialize()
      {
        if (parameter_error_count_ > 0 || verify() > 0)
        {
          Log::error() << "Reecb: cannot initialize with invalid configuration\n";
          return 1;
        }

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

        y_[VT] = std::sqrt(vr * vr + vi * vi);
        if (!Vref0_given_)
        {
          Vref0_ = static_cast<RealT>(y_[VT]);
        }

        y_[VMEAS]     = y_[VT];
        y_[VMEASSAFE] = Math::max(y_[VMEAS], static_cast<RealT>(0.01));

        const ScalarT ipcmd0 = y_[IPCMD];
        const ScalarT iqcmd0 = y_[IQCMD];
        const ScalarT pe0    = ipcmd0 * y_[VMEASSAFE];
        const ScalarT qgen0  = iqcmd0 * y_[VMEASSAFE];

        const ScalarT qref0 = qgen0;
        const ScalarT pref0 = Math::clamp(pe0, Pmin_, Pmax_);
        ScalarT       pfaref0{ZERO<RealT>};
        if (std::abs(static_cast<RealT>(pe0)) > ZERO<RealT>)
        {
          pfaref0 = static_cast<ScalarT>(std::atan(static_cast<RealT>(qref0 / pe0)));
        }

        pe_set_     = pe0;
        qgen_set_   = qgen0;
        qext_set_   = qref0;
        pfaref_set_ = pfaref0;
        pref_set_   = pref0;

        if (signals_.template isAttached<ReecbExternalVariables::QEXT>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::QEXT>(qref0);
        }
        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::PFAREF>(pfaref0);
        }
        if (signals_.template isAttached<ReecbExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::PREF>(pref0);
        }

        y_[PMEAS]              = pe0;
        const ScalarT qref_pf  = PfFlag_ * y_[PMEAS] * std::tan(pfaref_set_);
        const ScalarT qref_ext = pf_off_ * qext_set_;

        y_[SDIP]  = Math::outside(y_[VT], Vdip_, Vup_);
        y_[VERR]  = Math::deadband2(Vref0_ - y_[VMEAS], dbd1_, dbd2_);
        y_[IQV]   = Math::clamp(kqv_ * y_[VERR], Iql1_, Iqh1_);
        y_[QREF]  = qref_pf + qref_ext;
        y_[EQ]    = Math::clamp(y_[QREF], Qmin_, Qmax_) - qgen_set_;
        y_[QV]    = y_[QREF] / y_[VMEASSAFE];
        y_[PORD]  = pref_set_;
        y_[FPORD] = ZERO<RealT>;
        y_[RPORD] = ZERO<RealT>;

        static constexpr RealT kSmallTol  = static_cast<RealT>(1.0e-9);
        static constexpr RealT kSatMargin = static_cast<RealT>(0.1);

        auto steadyPiOutput = [&](const ScalarT interior_target,
                                  const ScalarT error,
                                  const ScalarT integral_rate,
                                  const ScalarT lower,
                                  const ScalarT upper) -> ScalarT
        {
          const RealT error_value = static_cast<RealT>(error);
          const RealT rate_value  = static_cast<RealT>(integral_rate);

          if (std::abs(error_value) <= kSmallTol || std::abs(rate_value) <= kSmallTol)
          {
            return interior_target;
          }

          if (rate_value > ZERO<RealT>)
          {
            return upper + static_cast<ScalarT>(kSatMargin);
          }

          return lower - static_cast<ScalarT>(kSatMargin);
        };

        const ScalarT xpiq_rate   = Kqi_ * y_[EQ];
        const ScalarT vpiq_target = VFlag_ * y_[VMEAS] + v_off_ * y_[QREF];
        const ScalarT vpiq_arg    = steadyPiOutput(vpiq_target,
                                                y_[EQ],
                                                xpiq_rate,
                                                static_cast<ScalarT>(Vmin_),
                                                static_cast<ScalarT>(Vmax_));
        y_[VPIQ]                  = Math::clamp(vpiq_arg, Vmin_, Vmax_);
        y_[EPIV]                  = VFlag_ * y_[VPIQ] + v_off_ * y_[QREF] - y_[VMEAS];
        y_[XPIQ]                  = vpiq_arg - Kqp_ * y_[EQ];

        const ScalarT iqbase_target = qgen0 / y_[VMEASSAFE];
        const ScalarT xpiv_rate     = Kvi_ * y_[EPIV];
        const ScalarT ip_star       = y_[PORD] / y_[VMEASSAFE];

        auto initializeReactiveBase = [&]()
        {
          const ScalarT piv_lower = -y_[IQMAX];
          const ScalarT piv_upper = y_[IQMAX];
          const ScalarT piv_arg   = steadyPiOutput(iqbase_target,
                                                 y_[EPIV],
                                                 xpiv_rate,
                                                 piv_lower,
                                                 piv_upper);
          y_[IQBASE]              = Math::clamp(piv_arg, -y_[IQMAX], y_[IQMAX]);
          y_[XPIV]                = piv_arg - Kvp_ * y_[EPIV];
        };

        auto initializeReactiveCommand = [&]()
        {
          y_[IQRAW] = QFlag_ * y_[IQBASE] + q_off_ * y_[QV] + y_[SDIP] * y_[IQV];
          y_[IQCMD] = Math::clamp(y_[IQRAW], -y_[IQMAX], y_[IQMAX]);
        };

        if (Pqflag_ == ZERO<RealT>)
        {
          y_[IQCIRC] = Imax_;
          y_[IQMAX]  = Imax_;
          initializeReactiveBase();
          initializeReactiveCommand();

          const ScalarT ip_radicand = Imax_ * Imax_ - y_[IQCMD] * y_[IQCMD];
          if (static_cast<RealT>(ip_radicand) < ZERO<RealT>)
          {
            Log::error() << "Reecb: initial active-current circle radicand is negative\n";
            return 1;
          }
          y_[IPCIRC] = std::sqrt(ip_radicand);
          y_[IPMAX]  = y_[IPCIRC];
          y_[IPCMD]  = Math::clamp(ip_star, ZERO<RealT>, y_[IPMAX]);
        }
        else
        {
          y_[IPCIRC] = Imax_;
          y_[IPMAX]  = Imax_;
          y_[IPCMD]  = Math::clamp(ip_star, ZERO<RealT>, y_[IPMAX]);

          const ScalarT iq_radicand = Imax_ * Imax_ - y_[IPCMD] * y_[IPCMD];
          if (static_cast<RealT>(iq_radicand) < ZERO<RealT>)
          {
            Log::error() << "Reecb: initial reactive-current circle radicand is negative\n";
            return 1;
          }
          y_[IQCIRC] = std::sqrt(iq_radicand);
          y_[IQMAX]  = y_[IQCIRC];
          initializeReactiveBase();
          initializeReactiveCommand();
        }

        std::fill(yp_.begin(), yp_.end(), ZERO<RealT>);

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Reecb<ScalarT, IdxT>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(ReecbInternalVariables::VMEAS)] = (Trv_ > ZERO<RealT>);
        tag_[static_cast<size_t>(ReecbInternalVariables::PMEAS)] = (Tp_ > ZERO<RealT>);
        tag_[static_cast<size_t>(ReecbInternalVariables::XPIQ)]  = true;
        tag_[static_cast<size_t>(ReecbInternalVariables::XPIV)]  = true;
        tag_[static_cast<size_t>(ReecbInternalVariables::QV)]    = true;
        tag_[static_cast<size_t>(ReecbInternalVariables::PORD)]  = true;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Reecb<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT* y,
          ScalarT* yp,
          ScalarT* wb,
          ScalarT* ws,
          ScalarT* f)
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

        const ScalarT pe     = ws[PE];
        const ScalarT qgen   = ws[QGEN];
        const ScalarT qext   = ws[QEXT];
        const ScalarT pfaref = ws[PFAREF];
        const ScalarT pref   = ws[PREF];

        const ScalarT not_sdip          = ONE<RealT> - y[SDIP];
        const ScalarT xpiq_arg          = Kqp_ * y[EQ] + y[XPIQ];
        const ScalarT xpiv_arg          = Kvp_ * y[EPIV] + y[XPIV];
        const ScalarT xpiq_limited_rate = Math::antiwindup(xpiq_arg, Kqi_ * y[EQ], Vmin_, Vmax_);
        const ScalarT xpiv_limited_rate = Math::antiwindup(xpiv_arg, Kvi_ * y[EPIV], -y[IQMAX], y[IQMAX]);
        const ScalarT pord_limited_rate = Math::antiwindup(y[PORD], y[RPORD], Pmin_, Pmax_);
        const ScalarT qref_pf           = PfFlag_ * y[PMEAS] * std::tan(pfaref);
        const ScalarT qref_ext          = pf_off_ * qext;
        const ScalarT epiv_target       = VFlag_ * y[VPIQ] + v_off_ * y[QREF];
        const ScalarT iqcirc_current_sq = Pqflag_ * y[IPCMD] * y[IPCMD];
        const ScalarT ipcirc_current_sq = (ONE<RealT> - Pqflag_) * y[IQCMD] * y[IQCMD];
        const ScalarT iqmax_target      = (ONE<RealT> - Pqflag_) * Imax_ + Pqflag_ * y[IQCIRC];
        const ScalarT ipmax_target      = Pqflag_ * Imax_ + (ONE<RealT> - Pqflag_) * y[IPCIRC];
        const ScalarT iqraw_target      = QFlag_ * y[IQBASE] + q_off_ * y[QV] + y[SDIP] * y[IQV];

        f[VMEAS] = -Trv_ * yp[VMEAS] - y[VMEAS] + y[VT];
        f[PMEAS] = -Tp_ * yp[PMEAS] - y[PMEAS] + pe;
        f[XPIQ]  = -yp[XPIQ] + not_sdip * xpiq_limited_rate;
        f[XPIV]  = -yp[XPIV] + not_sdip * xpiv_limited_rate;
        f[QV]    = -Tiq_ * yp[QV] - not_sdip * y[QV] + not_sdip * y[QREF] / y[VMEASSAFE];
        f[PORD]  = -yp[PORD] + not_sdip * pord_limited_rate;

        f[VT]        = -y[VT] * y[VT] + vr * vr + vi * vi;
        f[VMEASSAFE] = -y[VMEASSAFE] + Math::max(y[VMEAS], static_cast<RealT>(0.01));
        f[SDIP]      = -y[SDIP] + Math::outside(y[VT], Vdip_, Vup_);
        f[VERR]      = -y[VERR] + Math::deadband2(Vref0_ - y[VMEAS], dbd1_, dbd2_);
        f[IQV]       = -y[IQV] + Math::clamp(kqv_ * y[VERR], Iql1_, Iqh1_);
        f[QREF]      = -y[QREF] + qref_pf + qref_ext;
        f[EQ]        = -y[EQ] + Math::clamp(y[QREF], Qmin_, Qmax_) - qgen;
        f[VPIQ]      = -y[VPIQ] + Math::clamp(xpiq_arg, Vmin_, Vmax_);
        f[EPIV]      = -y[EPIV] + epiv_target - y[VMEAS];
        f[FPORD]     = -Tpord_ * y[FPORD] + pref - y[PORD];
        f[RPORD]     = -y[RPORD] + Math::clamp(y[FPORD], dPmin_, dPmax_);

        f[IQCIRC] = -y[IQCIRC] * y[IQCIRC] + Imax_ * Imax_ - iqcirc_current_sq;
        f[IPCIRC] = -y[IPCIRC] * y[IPCIRC] + Imax_ * Imax_ - ipcirc_current_sq;
        f[IQMAX]  = -y[IQMAX] + iqmax_target;
        f[IPMAX]  = -y[IPMAX] + ipmax_target;
        f[IQBASE] = -y[IQBASE] + Math::clamp(xpiv_arg, -y[IQMAX], y[IQMAX]);
        f[IQRAW]  = -y[IQRAW] + iqraw_target;
        f[IQCMD]  = -y[IQCMD] + Math::clamp(y[IQRAW], -y[IQMAX], y[IQMAX]);
        f[IPCMD]  = -y[IPCMD] + Math::clamp(y[PORD] / y[VMEASSAFE], ZERO<RealT>, y[IPMAX]);

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Reecb<ScalarT, IdxT>::evaluateResidual()
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
          ws_[PE] = toComponentBase(
              signals_.template readExternalVariable<ReecbExternalVariables::PE>());
          ws_indices_[PE] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PE>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QGEN>())
        {
          ws_[QGEN] = toComponentBase(
              signals_.template readExternalVariable<ReecbExternalVariables::QGEN>());
          ws_indices_[QGEN] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::QGEN>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QEXT>())
        {
          ws_[QEXT] = signals_.template readExternalVariable<ReecbExternalVariables::QEXT>();
          ws_indices_[QEXT] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::QEXT>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>())
        {
          ws_[PFAREF] = signals_.template readExternalVariable<ReecbExternalVariables::PFAREF>();
          ws_indices_[PFAREF] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PFAREF>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PREF>())
        {
          ws_[PREF] = signals_.template readExternalVariable<ReecbExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PREF>();
        }

        wb_[0] = Vr();
        wb_[1] = Vi();

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
        return 0;
      }
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

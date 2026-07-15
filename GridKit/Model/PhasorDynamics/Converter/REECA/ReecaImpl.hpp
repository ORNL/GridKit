/**
 * @file ReecaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REECA electrical-control model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECA/Reeca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECA/ReecaData.hpp>
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
      Reeca<scalar_type, index_type>::Reeca(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(ReecaInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      Reeca<scalar_type, index_type>::Reeca(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(ReecaInternalVariables::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Reeca<scalar_type, index_type>::~Reeca()
      {
      }

      template <typename scalar_type, typename index_type>
      scalar_type& Reeca<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      template <typename scalar_type, typename index_type>
      scalar_type& Reeca<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

      template <typename scalar_type, typename index_type>
      void Reeca<scalar_type, index_type>::setDerivedParameters()
      {
        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);
        Trv_eff_           = std::max(Trv_, TIME_CONSTANT_MINIMUM);
        Tp_eff_            = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        pf_off_            = ONE<RealT> - PfFlag_;
        v_off_             = ONE<RealT> - VFlag_;
        q_off_             = ONE<RealT> - QFlag_;
        p_off_             = ONE<RealT> - PFlag_;
        pq_off_            = ONE<RealT> - Pqflag_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Reeca<scalar_type, index_type>::toComponentBase(
          scalar_type value) const
      {
        return value * va_system_base_ / va_converter_base_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Reeca<scalar_type, index_type>::toSystemBase(
          scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      template <typename scalar_type, typename index_type>
      scalar_type Reeca<scalar_type, index_type>::currentLimit(
          scalar_type vdl_limit,
          scalar_type circle_limit) const
      {
        return vdl_limit * circle_limit / Math::max(vdl_limit, circle_limit);
      }

      template <typename scalar_type, typename index_type>
      void Reeca<scalar_type, index_type>::initModelParams(const ModelDataT& data)
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
            Log::error() << "Reeca: parameter '" << name << "' must be numeric\n";
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
          else if (const auto* real_value = std::get_if<RealT>(&value))
          {
            if (*real_value == ZERO<RealT> || *real_value == ONE<RealT>)
            {
              target = *real_value;
            }
            else
            {
              Log::error() << "Reeca: parameter '" << name << "' must be bool or 0/1\n";
              ++parameter_error_count_;
            }
          }
          else
          {
            Log::error() << "Reeca: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        if (!data.parameters.contains(Params::mva))
        {
          Log::error() << "Reeca: missing required parameter 'mva'\n";
          ++parameter_error_count_;
        }
        load_real(Params::mva, mva_base_, "mva");
        load_switch(Params::PfFlag, PfFlag_, "PfFlag");
        load_switch(Params::VFlag, VFlag_, "VFlag");
        load_switch(Params::QFlag, QFlag_, "QFlag");
        load_switch(Params::PFlag, PFlag_, "PFlag");
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
        load_real(Params::Iqfrz, Iqfrz_, "Iqfrz");
        load_real(Params::Thld, Thld_, "Thld");
        load_real(Params::Thld2, Thld2_, "Thld2");
        load_real(Params::Qmax, Qmax_, "Qmax");
        load_real(Params::Qmin, Qmin_, "Qmin");
        load_real(Params::Kqp, Kqp_, "Kqp");
        load_real(Params::Kqi, Kqi_, "Kqi");
        load_real(Params::Vmax, Vmax_, "Vmax");
        load_real(Params::Vmin, Vmin_, "Vmin");
        load_real(Params::Vref1, Vref1_, "Vref1");
        load_real(Params::Kvp, Kvp_, "Kvp");
        load_real(Params::Kvi, Kvi_, "Kvi");
        load_real(Params::Tiq, Tiq_, "Tiq");
        load_real(Params::Tpord, Tpord_, "Tpord");
        load_real(Params::dPmax, dPmax_, "dPmax");
        load_real(Params::dPmin, dPmin_, "dPmin");
        load_real(Params::Pmax, Pmax_, "Pmax");
        load_real(Params::Pmin, Pmin_, "Pmin");
        load_real(Params::Imax, Imax_, "Imax");
        load_real(Params::Vq1, Vq1_, "Vq1");
        load_real(Params::Iq1, Iq1_, "Iq1");
        load_real(Params::Vq2, Vq2_, "Vq2");
        load_real(Params::Iq2, Iq2_, "Iq2");
        load_real(Params::Vq3, Vq3_, "Vq3");
        load_real(Params::Iq3, Iq3_, "Iq3");
        load_real(Params::Vq4, Vq4_, "Vq4");
        load_real(Params::Iq4, Iq4_, "Iq4");
        load_real(Params::Vp1, Vp1_, "Vp1");
        load_real(Params::Ip1, Ip1_, "Ip1");
        load_real(Params::Vp2, Vp2_, "Vp2");
        load_real(Params::Ip2, Ip2_, "Ip2");
        load_real(Params::Vp3, Vp3_, "Vp3");
        load_real(Params::Ip3, Ip3_, "Ip3");
        load_real(Params::Vp4, Vp4_, "Vp4");
        load_real(Params::Ip4, Ip4_, "Ip4");
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Reeca<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Reeca<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        auto index     = [](ReecaInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::iqcmd, [this, index]
                      { return y_.getData()[index(ReecaInternalVariables::IQCMD)]; });
        monitor_->set(Variable::ipcmd, [this, index]
                      { return y_.getData()[index(ReecaInternalVariables::IPCMD)]; });
        monitor_->set(Variable::vmeas, [this, index]
                      { return y_.getData()[index(ReecaInternalVariables::VMEAS)]; });
        monitor_->set(Variable::pmeas, [this, index]
                      { return y_.getData()[index(ReecaInternalVariables::PMEAS)]; });
      }

      template <typename scalar_type, typename index_type>
      int Reeca<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reeca<scalar_type, index_type>::allocate()
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

        auto signal_size = static_cast<size_t>(ReecaExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<ReecaInternalVariables::IQCMD>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<ReecaInternalVariables::IQCMD>()->set(
              &y[static_cast<size_t>(ReecaInternalVariables::IQCMD)],
              &(this->getVariableIndex(static_cast<IdxT>(ReecaInternalVariables::IQCMD))));
        }

        if (signals_.template isAssigned<ReecaInternalVariables::IPCMD>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<ReecaInternalVariables::IPCMD>()->set(
              &y[static_cast<size_t>(ReecaInternalVariables::IPCMD)],
              &(this->getVariableIndex(static_cast<IdxT>(ReecaInternalVariables::IPCMD))));
        }

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reeca<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Reeca: " << message << '\n';
            ret += 1;
          }
        };

        if (bus_ == nullptr)
        {
          Log::error() << "Reeca: bus pointer is null\n";
          ret += 1;
        }

        check(mva_base_ > ZERO<RealT>, "mva must be positive");
        check(va_converter_base_ > ZERO<RealT>, "converter VA base must be positive");
        check(PfFlag_ == ZERO<RealT> || PfFlag_ == ONE<RealT>, "PfFlag must be 0 or 1");
        check(VFlag_ == ZERO<RealT> || VFlag_ == ONE<RealT>, "VFlag must be 0 or 1");
        check(QFlag_ == ZERO<RealT> || QFlag_ == ONE<RealT>, "QFlag must be 0 or 1");
        check(PFlag_ == ZERO<RealT> || PFlag_ == ONE<RealT>, "PFlag must be 0 or 1");
        check(Pqflag_ == ZERO<RealT> || Pqflag_ == ONE<RealT>, "Pqflag must be 0 or 1");
        check(Trv_ >= ZERO<RealT>, "Trv must be non-negative");
        check(Tp_ >= ZERO<RealT>, "Tp must be non-negative");
        check(ZERO<RealT> <= Vdip_ && Vdip_ < Vup_, "Vdip/Vup must satisfy 0 <= Vdip < Vup");
        check(dbd1_ <= ZERO<RealT> && ZERO<RealT> <= dbd2_, "dbd1 <= 0 <= dbd2 is required");
        check(Iql1_ <= Iqh1_, "Iql1 must be less than or equal to Iqh1");
        check(Thld_ == ZERO<RealT>, "Thld must be zero in this REECA implementation");
        check(Thld2_ == ZERO<RealT>, "Thld2 must be zero in this REECA implementation");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax");
        check(Tiq_ > ZERO<RealT>, "Tiq must be positive");
        check(Tpord_ > ZERO<RealT>, "Tpord must be positive");
        check(dPmin_ < ZERO<RealT> && ZERO<RealT> < dPmax_, "dPmin < 0 < dPmax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");
        check(Imax_ >= ZERO<RealT>, "Imax must be non-negative");
        check(ZERO<RealT> <= Vq1_ && Vq1_ < Vq2_ && Vq2_ < Vq3_ && Vq3_ < Vq4_,
              "Vq breakpoints must satisfy 0 <= Vq1 < Vq2 < Vq3 < Vq4");
        check(Iq1_ >= ZERO<RealT> && Iq2_ >= ZERO<RealT> && Iq3_ >= ZERO<RealT>
                  && Iq4_ >= ZERO<RealT>,
              "Iq VDL limits must be non-negative");
        check(ZERO<RealT> <= Vp1_ && Vp1_ < Vp2_ && Vp2_ < Vp3_ && Vp3_ < Vp4_,
              "Vp breakpoints must satisfy 0 <= Vp1 < Vp2 < Vp3 < Vp4");
        check(Ip1_ >= ZERO<RealT> && Ip2_ >= ZERO<RealT> && Ip3_ >= ZERO<RealT>
                  && Ip4_ >= ZERO<RealT>,
              "Ip VDL limits must be non-negative");
        check(PFlag_ != ONE<RealT>
                  || signals_.template isAttached<ReecaExternalVariables::OMEGA>(),
              "omegag is required when PFlag is 1");

        const bool has_pe   = signals_.template isAttached<ReecaExternalVariables::PE>();
        const bool has_qgen = signals_.template isAttached<ReecaExternalVariables::QGEN>();

        if (has_pe != has_qgen)
        {
          Log::error() << "Reeca: pe and qgen must be connected together\n";
          ret += 1;
        }

        auto check_attached_signal = [&](bool attached, bool linked, const char* name)
        {
          if (attached && !linked)
          {
            Log::error() << "Reeca: " << name << " signal attached with no linked variable\n";
            ret += 1;
          }
        };

        check_attached_signal(
            signals_.template isAttached<ReecaExternalVariables::PE>(),
            signals_.template isAttached<ReecaExternalVariables::PE>()
                && signals_.template isLinked<ReecaExternalVariables::PE>(),
            "pe");
        check_attached_signal(
            signals_.template isAttached<ReecaExternalVariables::QGEN>(),
            signals_.template isAttached<ReecaExternalVariables::QGEN>()
                && signals_.template isLinked<ReecaExternalVariables::QGEN>(),
            "qgen");
        check_attached_signal(
            signals_.template isAttached<ReecaExternalVariables::OMEGA>(),
            signals_.template isAttached<ReecaExternalVariables::OMEGA>()
                && signals_.template isLinked<ReecaExternalVariables::OMEGA>(),
            "omegag");
        check_attached_signal(
            signals_.template isAttached<ReecaExternalVariables::QEXT>(),
            signals_.template isAttached<ReecaExternalVariables::QEXT>()
                && signals_.template isLinked<ReecaExternalVariables::QEXT>(),
            "qext");
        check_attached_signal(
            signals_.template isAttached<ReecaExternalVariables::PFAREF>(),
            signals_.template isAttached<ReecaExternalVariables::PFAREF>()
                && signals_.template isLinked<ReecaExternalVariables::PFAREF>(),
            "pfaref");
        check_attached_signal(
            signals_.template isAttached<ReecaExternalVariables::PREF>(),
            signals_.template isAttached<ReecaExternalVariables::PREF>()
                && signals_.template isLinked<ReecaExternalVariables::PREF>(),
            "pref");

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int Reeca<scalar_type, index_type>::initialize()
      {
        if (parameter_error_count_ > 0 || verify() > 0)
        {
          Log::error() << "Reeca: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        const auto VMEAS     = static_cast<size_t>(ReecaInternalVariables::VMEAS);
        const auto PMEAS     = static_cast<size_t>(ReecaInternalVariables::PMEAS);
        const auto XPIQ      = static_cast<size_t>(ReecaInternalVariables::XPIQ);
        const auto XPIV      = static_cast<size_t>(ReecaInternalVariables::XPIV);
        const auto QV        = static_cast<size_t>(ReecaInternalVariables::QV);
        const auto PORD      = static_cast<size_t>(ReecaInternalVariables::PORD);
        const auto VT        = static_cast<size_t>(ReecaInternalVariables::VT);
        const auto VMEASSAFE = static_cast<size_t>(ReecaInternalVariables::VMEASSAFE);
        const auto SDIP      = static_cast<size_t>(ReecaInternalVariables::SDIP);
        const auto VERR      = static_cast<size_t>(ReecaInternalVariables::VERR);
        const auto IQV       = static_cast<size_t>(ReecaInternalVariables::IQV);
        const auto QREF      = static_cast<size_t>(ReecaInternalVariables::QREF);
        const auto EQ        = static_cast<size_t>(ReecaInternalVariables::EQ);
        const auto VPIQ      = static_cast<size_t>(ReecaInternalVariables::VPIQ);
        const auto EPIV      = static_cast<size_t>(ReecaInternalVariables::EPIV);
        const auto FPORD     = static_cast<size_t>(ReecaInternalVariables::FPORD);
        const auto RPORD     = static_cast<size_t>(ReecaInternalVariables::RPORD);
        const auto IQCIRC    = static_cast<size_t>(ReecaInternalVariables::IQCIRC);
        const auto IPCIRC    = static_cast<size_t>(ReecaInternalVariables::IPCIRC);
        const auto IQMAX     = static_cast<size_t>(ReecaInternalVariables::IQMAX);
        const auto IPMAX     = static_cast<size_t>(ReecaInternalVariables::IPMAX);
        const auto IQBASE    = static_cast<size_t>(ReecaInternalVariables::IQBASE);
        const auto IQRAW     = static_cast<size_t>(ReecaInternalVariables::IQRAW);
        const auto IQCMD     = static_cast<size_t>(ReecaInternalVariables::IQCMD);
        const auto IPCMD     = static_cast<size_t>(ReecaInternalVariables::IPCMD);

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();

        y[VT] = std::sqrt(vr * vr + vi * vi);

        const RealT vt0 = static_cast<RealT>(y[VT]);
        if (vt0 <= ZERO<RealT>)
        {
          Log::error() << "Reeca: terminal voltage magnitude must be positive at initialization\n";
          return 1;
        }

        if (!(Vdip_ < vt0 && vt0 < Vup_))
        {
          Log::error() << "Reeca: standard initialization requires Vdip < Vt0 < Vup\n";
          return 1;
        }

        if (!Vref0_given_)
        {
          Vref0_ = vt0;
        }

        y[VMEAS]     = y[VT];
        y[VMEASSAFE] = Math::max(y[VMEAS], VMEAS_MINIMUM);

        const ScalarT ipcmd0 = toComponentBase(y[IPCMD]);
        const ScalarT iqcmd0 = toComponentBase(y[IQCMD]);
        const ScalarT pe0    = ipcmd0 * y[VMEASSAFE];
        const ScalarT qgen0  = iqcmd0 * y[VMEASSAFE];

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<ReecaExternalVariables::OMEGA>())
        {
          omega0 = signals_.template readExternalVariable<ReecaExternalVariables::OMEGA>();
        }

        const ScalarT pref_denominator = p_off_ + PFlag_ * omega0;
        if (static_cast<RealT>(pref_denominator) == ZERO<RealT>)
        {
          Log::error() << "Reeca: pref initialization denominator is zero\n";
          return 1;
        }

        const ScalarT qext0 = qgen0;
        const ScalarT pref0 = pe0 / pref_denominator;
        ScalarT       pfaref0{ZERO<RealT>};
        if (std::abs(static_cast<RealT>(pe0)) > INIT_TOL)
        {
          pfaref0 = static_cast<ScalarT>(std::atan(static_cast<RealT>(qgen0 / pe0)));
        }

        pe_set_     = toSystemBase(pe0);
        qgen_set_   = toSystemBase(qgen0);
        omega_set_  = omega0;
        qext_set_   = toSystemBase(qext0);
        pfaref_set_ = pfaref0;
        pref_set_   = toSystemBase(pref0);

        if (signals_.template isAttached<ReecaExternalVariables::PE>())
        {
          signals_.template writeExternalVariable<ReecaExternalVariables::PE>(pe_set_);
        }
        if (signals_.template isAttached<ReecaExternalVariables::QGEN>())
        {
          signals_.template writeExternalVariable<ReecaExternalVariables::QGEN>(qgen_set_);
        }
        if (signals_.template isAttached<ReecaExternalVariables::QEXT>())
        {
          signals_.template writeExternalVariable<ReecaExternalVariables::QEXT>(qext_set_);
        }
        if (signals_.template isAttached<ReecaExternalVariables::PFAREF>())
        {
          signals_.template writeExternalVariable<ReecaExternalVariables::PFAREF>(pfaref0);
        }
        if (signals_.template isAttached<ReecaExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<ReecaExternalVariables::PREF>(pref_set_);
        }

        y[PMEAS]               = pe0;
        const ScalarT qref_pf  = PfFlag_ * y[PMEAS] * std::tan(pfaref_set_);
        const ScalarT qref_ext = pf_off_ * qext0;

        y[SDIP] = Math::inside(y[VT], Vdip_, Vup_);
        y[VERR] = Math::deadband2(Vref0_ - y[VMEAS], dbd1_, dbd2_);
        y[IQV]  = Math::clamp(kqv_ * y[VERR], Iql1_, Iqh1_);
        y[QREF] = qref_pf + qref_ext;

        const RealT qref0_value = static_cast<RealT>(y[QREF]);
        if (qref0_value < Qmin_ || qref0_value > Qmax_)
        {
          Log::error() << "Reeca: standard initialization requires Qref0 within Qmin/Qmax\n";
          return 1;
        }

        y[EQ] = Math::clamp(y[QREF], Qmin_, Qmax_) - qgen0;
        y[QV] = y[QREF] / y[VMEASSAFE];

        const ScalarT pord0       = pref_denominator * pref0;
        const RealT   pord0_value = static_cast<RealT>(pord0);
        if (pord0_value < Pmin_ || pord0_value > Pmax_)
        {
          Log::error() << "Reeca: standard initialization requires Pord0 within Pmin/Pmax\n";
          return 1;
        }

        y[PORD]  = pord0;
        y[FPORD] = ZERO<RealT>;
        y[RPORD] = ZERO<RealT>;

        auto awinit = [](const ScalarT target,
                         const ScalarT rate,
                         const ScalarT lower,
                         const ScalarT upper) -> ScalarT
        {
          if (std::abs(static_cast<RealT>(rate)) <= INIT_TOL)
          {
            return target;
          }
          return rate > ZERO<RealT> ? upper + static_cast<ScalarT>(SAT_MARGIN)
                                    : lower - static_cast<ScalarT>(SAT_MARGIN);
        };

        const ScalarT xpiq_rate   = Kqi_ * y[EQ];
        const ScalarT vpiq_target = VFlag_ * y[VMEAS] + v_off_ * (y[QREF] + Vref1_);
        const ScalarT vpiq_arg    = awinit(vpiq_target,
                                        xpiq_rate,
                                        static_cast<ScalarT>(Vmin_),
                                        static_cast<ScalarT>(Vmax_));
        y[VPIQ]                   = Math::clamp(vpiq_arg, Vmin_, Vmax_);
        y[EPIV]                   = VFlag_ * y[VPIQ] + v_off_ * (y[QREF] + Vref1_) - y[VMEAS];
        y[XPIQ]                   = vpiq_arg - Kqp_ * y[EQ];

        const ScalarT gq0 = static_cast<ScalarT>(Iq1_)
                            + Math::linseg(y[VMEAS], Vq1_, Vq2_, Iq2_ - Iq1_)
                            + Math::linseg(y[VMEAS], Vq2_, Vq3_, Iq3_ - Iq2_)
                            + Math::linseg(y[VMEAS], Vq3_, Vq4_, Iq4_ - Iq3_);
        const ScalarT gp0 = static_cast<ScalarT>(Ip1_)
                            + Math::linseg(y[VMEAS], Vp1_, Vp2_, Ip2_ - Ip1_)
                            + Math::linseg(y[VMEAS], Vp2_, Vp3_, Ip3_ - Ip2_)
                            + Math::linseg(y[VMEAS], Vp3_, Vp4_, Ip4_ - Ip3_);

        const ScalarT iqbase_target = qgen0 / y[VMEASSAFE];
        const ScalarT xpiv_rate     = Kvi_ * y[EPIV];
        const ScalarT ip_star       = y[PORD] / y[VMEASSAFE];

        auto initializeReactiveBase = [&]()
        {
          const ScalarT piv_arg = awinit(iqbase_target, xpiv_rate, -y[IQMAX], y[IQMAX]);
          y[IQBASE]             = Math::clamp(piv_arg, -y[IQMAX], y[IQMAX]);
          y[XPIV]               = piv_arg - Kvp_ * y[EPIV];
        };

        auto initializeReactiveCommand = [&]()
        {
          y[IQRAW] = QFlag_ * y[IQBASE] + q_off_ * y[QV] + (ONE<RealT> - y[SDIP]) * y[IQV];
          y[IQCMD] = toSystemBase(Math::clamp(y[IQRAW], -y[IQMAX], y[IQMAX]));
        };

        auto rejectOutsideReactiveLimits = [&](const ScalarT value, const char* message)
        {
          const RealT limit = static_cast<RealT>(y[IQMAX]);
          if (static_cast<RealT>(value) < -limit || static_cast<RealT>(value) > limit)
          {
            Log::error() << "Reeca: " << message << '\n';
            return true;
          }
          return false;
        };

        auto rejectOutsideActiveLimits = [&](const ScalarT value, const char* message)
        {
          const RealT limit = static_cast<RealT>(y[IPMAX]);
          if (static_cast<RealT>(value) < ZERO<RealT> || static_cast<RealT>(value) > limit)
          {
            Log::error() << "Reeca: " << message << '\n';
            return true;
          }
          return false;
        };

        if (Pqflag_ == ZERO<RealT>)
        {
          y[IQCIRC] = Imax_;
          y[IQMAX]  = currentLimit(gq0, y[IQCIRC]);
          if (rejectOutsideReactiveLimits(
                  iqbase_target,
                  "standard initialization requires Iq0 within reactive-current limits"))
          {
            return 1;
          }
          initializeReactiveBase();
          initializeReactiveCommand();
          if (rejectOutsideReactiveLimits(
                  y[IQRAW],
                  "standard initialization requires raw Iq0 within reactive-current limits"))
          {
            return 1;
          }

          const ScalarT iqcmd       = toComponentBase(y[IQCMD]);
          const ScalarT ip_radicand = Imax_ * Imax_ - iqcmd * iqcmd;
          if (static_cast<RealT>(ip_radicand) < ZERO<RealT>)
          {
            Log::error() << "Reeca: initial active-current circle radicand is negative\n";
            return 1;
          }
          y[IPCIRC] = std::sqrt(ip_radicand);
          y[IPMAX]  = currentLimit(gp0, y[IPCIRC]);
          if (rejectOutsideActiveLimits(
                  ip_star,
                  "standard initialization requires Ip0 within active-current limits"))
          {
            return 1;
          }
          y[IPCMD] = toSystemBase(Math::clamp(ip_star, ZERO<RealT>, y[IPMAX]));
        }
        else
        {
          y[IPCIRC] = Imax_;
          y[IPMAX]  = currentLimit(gp0, y[IPCIRC]);
          if (rejectOutsideActiveLimits(
                  ip_star,
                  "standard initialization requires Ip0 within active-current limits"))
          {
            return 1;
          }
          y[IPCMD] = toSystemBase(Math::clamp(ip_star, ZERO<RealT>, y[IPMAX]));

          const ScalarT ipcmd       = toComponentBase(y[IPCMD]);
          const ScalarT iq_radicand = Imax_ * Imax_ - ipcmd * ipcmd;
          if (static_cast<RealT>(iq_radicand) < ZERO<RealT>)
          {
            Log::error() << "Reeca: initial reactive-current circle radicand is negative\n";
            return 1;
          }
          y[IQCIRC] = std::sqrt(iq_radicand);
          y[IQMAX]  = currentLimit(gq0, y[IQCIRC]);
          if (rejectOutsideReactiveLimits(
                  iqbase_target,
                  "standard initialization requires Iq0 within reactive-current limits"))
          {
            return 1;
          }
          initializeReactiveBase();
          initializeReactiveCommand();
          if (rejectOutsideReactiveLimits(
                  y[IQRAW],
                  "standard initialization requires raw Iq0 within reactive-current limits"))
          {
            return 1;
          }
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
      int Reeca<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(ReecaInternalVariables::VMEAS)] = true;
        tag_[static_cast<size_t>(ReecaInternalVariables::PMEAS)] = true;
        tag_[static_cast<size_t>(ReecaInternalVariables::XPIQ)]  = true;
        tag_[static_cast<size_t>(ReecaInternalVariables::XPIV)]  = true;
        tag_[static_cast<size_t>(ReecaInternalVariables::QV)]    = true;
        tag_[static_cast<size_t>(ReecaInternalVariables::PORD)]  = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reeca<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Reeca<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto VMEAS     = static_cast<size_t>(ReecaInternalVariables::VMEAS);
        const auto PMEAS     = static_cast<size_t>(ReecaInternalVariables::PMEAS);
        const auto XPIQ      = static_cast<size_t>(ReecaInternalVariables::XPIQ);
        const auto XPIV      = static_cast<size_t>(ReecaInternalVariables::XPIV);
        const auto QV        = static_cast<size_t>(ReecaInternalVariables::QV);
        const auto PORD      = static_cast<size_t>(ReecaInternalVariables::PORD);
        const auto VT        = static_cast<size_t>(ReecaInternalVariables::VT);
        const auto VMEASSAFE = static_cast<size_t>(ReecaInternalVariables::VMEASSAFE);
        const auto SDIP      = static_cast<size_t>(ReecaInternalVariables::SDIP);
        const auto VERR      = static_cast<size_t>(ReecaInternalVariables::VERR);
        const auto IQV       = static_cast<size_t>(ReecaInternalVariables::IQV);
        const auto QREF      = static_cast<size_t>(ReecaInternalVariables::QREF);
        const auto EQ        = static_cast<size_t>(ReecaInternalVariables::EQ);
        const auto VPIQ      = static_cast<size_t>(ReecaInternalVariables::VPIQ);
        const auto EPIV      = static_cast<size_t>(ReecaInternalVariables::EPIV);
        const auto FPORD     = static_cast<size_t>(ReecaInternalVariables::FPORD);
        const auto RPORD     = static_cast<size_t>(ReecaInternalVariables::RPORD);
        const auto IQCIRC    = static_cast<size_t>(ReecaInternalVariables::IQCIRC);
        const auto IPCIRC    = static_cast<size_t>(ReecaInternalVariables::IPCIRC);
        const auto IQMAX     = static_cast<size_t>(ReecaInternalVariables::IQMAX);
        const auto IPMAX     = static_cast<size_t>(ReecaInternalVariables::IPMAX);
        const auto IQBASE    = static_cast<size_t>(ReecaInternalVariables::IQBASE);
        const auto IQRAW     = static_cast<size_t>(ReecaInternalVariables::IQRAW);
        const auto IQCMD     = static_cast<size_t>(ReecaInternalVariables::IQCMD);
        const auto IPCMD     = static_cast<size_t>(ReecaInternalVariables::IPCMD);

        const auto PE     = static_cast<size_t>(ReecaExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecaExternalVariables::QGEN);
        const auto OMEGA  = static_cast<size_t>(ReecaExternalVariables::OMEGA);
        const auto QEXT   = static_cast<size_t>(ReecaExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecaExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecaExternalVariables::PREF);

        const ScalarT vmeas      = y[VMEAS];
        const ScalarT pmeas      = y[PMEAS];
        const ScalarT xpiq       = y[XPIQ];
        const ScalarT xpiv       = y[XPIV];
        const ScalarT qv         = y[QV];
        const ScalarT pord       = y[PORD];
        const ScalarT vt         = y[VT];
        const ScalarT vmeas_safe = y[VMEASSAFE];
        const ScalarT sdip       = y[SDIP];
        const ScalarT verr       = y[VERR];
        const ScalarT iqv        = y[IQV];
        const ScalarT qref       = y[QREF];
        const ScalarT eq         = y[EQ];
        const ScalarT vpiq       = y[VPIQ];
        const ScalarT epiv       = y[EPIV];
        const ScalarT fpord      = y[FPORD];
        const ScalarT rpord      = y[RPORD];
        const ScalarT iqcirc     = y[IQCIRC];
        const ScalarT ipcirc     = y[IPCIRC];
        const ScalarT iqmax      = y[IQMAX];
        const ScalarT ipmax      = y[IPMAX];
        const ScalarT iqbase     = y[IQBASE];
        const ScalarT iqraw      = y[IQRAW];
        const ScalarT iqcmd      = y[IQCMD];
        const ScalarT ipcmd      = y[IPCMD];

        const ScalarT vmeas_dot = yp[VMEAS];
        const ScalarT pmeas_dot = yp[PMEAS];
        const ScalarT xpiq_dot  = yp[XPIQ];
        const ScalarT xpiv_dot  = yp[XPIV];
        const ScalarT qv_dot    = yp[QV];
        const ScalarT pord_dot  = yp[PORD];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT pe              = toComponentBase(ws[PE]);
        const ScalarT qgen            = toComponentBase(ws[QGEN]);
        const ScalarT omega           = ws[OMEGA];
        const ScalarT qext            = toComponentBase(ws[QEXT]);
        const ScalarT pfaref          = ws[PFAREF];
        const ScalarT pref            = toComponentBase(ws[PREF]);
        const ScalarT iqcmd_component = toComponentBase(iqcmd);
        const ScalarT ipcmd_component = toComponentBase(ipcmd);

        const ScalarT xpiq_arg          = Kqp_ * eq + xpiq;
        const ScalarT xpiv_arg          = Kvp_ * epiv + xpiv;
        const ScalarT xpiq_limited_rate = Math::antiwindup(xpiq_arg, Kqi_ * eq, Vmin_, Vmax_);
        const ScalarT xpiv_limited_rate =
            Math::antiwindup(xpiv_arg, Kvi_ * epiv, -iqmax, iqmax);
        const ScalarT pord_limited_rate = Math::antiwindup(pord, rpord, Pmin_, Pmax_);
        const ScalarT qref_pf           = PfFlag_ * pmeas * std::tan(pfaref);
        const ScalarT qref_ext          = pf_off_ * qext;
        const ScalarT epiv_target       = VFlag_ * vpiq + v_off_ * (qref + Vref1_);
        const ScalarT iqcirc_current_sq = Pqflag_ * ipcmd_component * ipcmd_component;
        const ScalarT ipcirc_current_sq = pq_off_ * iqcmd_component * iqcmd_component;
        const ScalarT gq                = static_cast<ScalarT>(Iq1_)
                           + Math::linseg(vmeas, Vq1_, Vq2_, Iq2_ - Iq1_)
                           + Math::linseg(vmeas, Vq2_, Vq3_, Iq3_ - Iq2_)
                           + Math::linseg(vmeas, Vq3_, Vq4_, Iq4_ - Iq3_);
        const ScalarT gp = static_cast<ScalarT>(Ip1_)
                           + Math::linseg(vmeas, Vp1_, Vp2_, Ip2_ - Ip1_)
                           + Math::linseg(vmeas, Vp2_, Vp3_, Ip3_ - Ip2_)
                           + Math::linseg(vmeas, Vp3_, Vp4_, Ip4_ - Ip3_);
        const ScalarT iqraw_target =
            QFlag_ * iqbase + q_off_ * qv + (ONE<RealT> - sdip) * iqv;
        const ScalarT iq_limit = currentLimit(gq, iqcirc);
        const ScalarT ip_limit = currentLimit(gp, ipcirc);

        f[VMEAS] = -vmeas_dot + (vt - vmeas) / Trv_eff_;
        f[PMEAS] = -pmeas_dot + (pe - pmeas) / Tp_eff_;
        f[XPIQ]  = -xpiq_dot + sdip * xpiq_limited_rate;
        f[XPIV]  = -xpiv_dot + sdip * xpiv_limited_rate;
        f[QV]    = -qv_dot + sdip * (qref / vmeas_safe - qv) / Tiq_;
        f[PORD]  = -pord_dot + sdip * pord_limited_rate;

        f[VT]        = -vt * vt + vr * vr + vi * vi;
        f[VMEASSAFE] = -vmeas_safe + Math::max(vmeas, VMEAS_MINIMUM);
        f[SDIP]      = -sdip + Math::inside(vt, Vdip_, Vup_);
        f[VERR]      = -verr + Math::deadband2(Vref0_ - vmeas, dbd1_, dbd2_);
        f[IQV]       = -iqv + Math::clamp(kqv_ * verr, Iql1_, Iqh1_);
        f[QREF]      = -qref + qref_pf + qref_ext;
        f[EQ]        = -eq + Math::clamp(qref, Qmin_, Qmax_) - qgen;
        f[VPIQ]      = -vpiq + Math::clamp(xpiq_arg, Vmin_, Vmax_);
        f[EPIV]      = -epiv + epiv_target - vmeas;
        f[FPORD]     = -fpord + ((p_off_ + PFlag_ * omega) * pref - pord) / Tpord_;
        f[RPORD]     = -rpord + Math::clamp(fpord, dPmin_, dPmax_);

        f[IQCIRC] = -iqcirc * iqcirc + Imax_ * Imax_ - iqcirc_current_sq;
        f[IPCIRC] = -ipcirc * ipcirc + Imax_ * Imax_ - ipcirc_current_sq;
        f[IQMAX]  = -iqmax + iq_limit;
        f[IPMAX]  = -ipmax + ip_limit;
        f[IQBASE] = -iqbase + Math::clamp(xpiv_arg, -iqmax, iqmax);
        f[IQRAW]  = -iqraw + iqraw_target;
        f[IQCMD]  = -iqcmd + toSystemBase(Math::clamp(iqraw, -iqmax, iqmax));
        f[IPCMD]  = -ipcmd + toSystemBase(Math::clamp(pord / vmeas_safe, ZERO<RealT>, ipmax));

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reeca<scalar_type, index_type>::evaluateResidual()
      {
        const auto PE     = static_cast<size_t>(ReecaExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecaExternalVariables::QGEN);
        const auto OMEGA  = static_cast<size_t>(ReecaExternalVariables::OMEGA);
        const auto QEXT   = static_cast<size_t>(ReecaExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecaExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecaExternalVariables::PREF);

        ws_[PE]     = pe_set_;
        ws_[QGEN]   = qgen_set_;
        ws_[OMEGA]  = omega_set_;
        ws_[QEXT]   = qext_set_;
        ws_[PFAREF] = pfaref_set_;
        ws_[PREF]   = pref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<ReecaExternalVariables::PE>())
        {
          ws_[PE] = signals_.template readExternalVariable<ReecaExternalVariables::PE>();
          ws_indices_[PE] =
              signals_.template readExternalVariableIndex<ReecaExternalVariables::PE>();
        }
        if (signals_.template isAttached<ReecaExternalVariables::QGEN>())
        {
          ws_[QGEN] = signals_.template readExternalVariable<ReecaExternalVariables::QGEN>();
          ws_indices_[QGEN] =
              signals_.template readExternalVariableIndex<ReecaExternalVariables::QGEN>();
        }
        if (signals_.template isAttached<ReecaExternalVariables::OMEGA>())
        {
          ws_[OMEGA] = signals_.template readExternalVariable<ReecaExternalVariables::OMEGA>();
          ws_indices_[OMEGA] =
              signals_.template readExternalVariableIndex<ReecaExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<ReecaExternalVariables::QEXT>())
        {
          ws_[QEXT] = signals_.template readExternalVariable<ReecaExternalVariables::QEXT>();
          ws_indices_[QEXT] =
              signals_.template readExternalVariableIndex<ReecaExternalVariables::QEXT>();
        }
        if (signals_.template isAttached<ReecaExternalVariables::PFAREF>())
        {
          ws_[PFAREF] = signals_.template readExternalVariable<ReecaExternalVariables::PFAREF>();
          ws_indices_[PFAREF] =
              signals_.template readExternalVariableIndex<ReecaExternalVariables::PFAREF>();
        }
        if (signals_.template isAttached<ReecaExternalVariables::PREF>())
        {
          ws_[PREF] = signals_.template readExternalVariable<ReecaExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<ReecaExternalVariables::PREF>();
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

/**
 * @file RepcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REPCA phasor-dynamics plant-control model.
 */

#pragma once

#include <algorithm>
#include <cmath>
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
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(RepcaInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(RepcaInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::~Repca()
      {
      }

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
            Log::error() << "Repca: parameter '" << name << "' must be numeric\n";
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
            Log::error() << "Repca: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        load_real(Params::mva, mva_, "mva");
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
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::setDerivedParameters()
      {
        Tfltr_ = std::max(Tfltr_, TIME_CONSTANT_MINIMUM);
        Tp_    = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        Tlag_  = std::max(Tlag_, TIME_CONSTANT_MINIMUM);

        va_component_base_ = mva_ * static_cast<RealT>(1.0e6);

        vcomp_on_  = VcompFlag_ ? ONE<RealT> : ZERO<RealT>;
        vcomp_off_ = ONE<RealT> - vcomp_on_;
        ref_on_    = RefFlag_ ? ONE<RealT> : ZERO<RealT>;
        ref_off_   = ONE<RealT> - ref_on_;
        freq_on_   = Freqflag_ ? ONE<RealT> : ZERO<RealT>;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Repca<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * va_system_base_ / va_component_base_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Repca<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Repca<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        auto index     = [](RepcaInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::qext, [this, index]
                      { return y_.getData()[index(RepcaInternalVariables::QEXT)]; });
        monitor_->set(Variable::pext, [this, index]
                      { return y_.getData()[index(RepcaInternalVariables::PEXT)]; });
        monitor_->set(Variable::vmeas, [this, index]
                      { return y_.getData()[index(RepcaInternalVariables::VMEAS)]; });
        monitor_->set(Variable::qmeas, [this, index]
                      { return y_.getData()[index(RepcaInternalVariables::QMEAS)]; });
        monitor_->set(Variable::pmeas, [this, index]
                      { return y_.getData()[index(RepcaInternalVariables::PMEAS)]; });
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::allocate()
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

        auto signal_size = static_cast<size_t>(RepcaExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<RepcaInternalVariables::QEXT>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<RepcaInternalVariables::QEXT>()->set(
              &y[static_cast<size_t>(RepcaInternalVariables::QEXT)],
              &(this->getVariableIndex(static_cast<IdxT>(RepcaInternalVariables::QEXT))));
        }

        if (signals_.template isAssigned<RepcaInternalVariables::PEXT>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<RepcaInternalVariables::PEXT>()->set(
              &y[static_cast<size_t>(RepcaInternalVariables::PEXT)],
              &(this->getVariableIndex(static_cast<IdxT>(RepcaInternalVariables::PEXT))));
        }

        allocated_ = true;
        return 0;
      }

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

        if (bus_ == nullptr)
        {
          Log::error() << "Repca: bus pointer is null\n";
          ret += 1;
        }

        check(mva_ > ZERO<RealT>, "mva must be positive");
        check(Tfv_ > ZERO<RealT>, "Tfv must be positive");
        check(dbdlow_ <= ZERO<RealT> && ZERO<RealT> <= dbdupper_,
              "dbdlow <= 0 <= dbdupper is required");
        check(emin_ <= ZERO<RealT> && ZERO<RealT> <= emax_,
              "emin <= 0 <= emax is required");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(fdbd1_ <= ZERO<RealT> && ZERO<RealT> <= fdbd2_,
              "fdbd1 <= 0 <= fdbd2 is required");
        check(femin_ <= ZERO<RealT> && ZERO<RealT> <= femax_,
              "femin <= 0 <= femax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");

        auto check_required_signal = [&](bool attached, bool linked, const char* name)
        {
          if (!attached)
          {
            Log::error() << "Repca: " << name << " signal is required\n";
            ret += 1;
          }
          else if (!linked)
          {
            Log::error() << "Repca: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_required_signal(
            signals_.template isAttached<RepcaExternalVariables::IBRANCHR>(),
            signals_.template isAttached<RepcaExternalVariables::IBRANCHR>()
                && signals_.template isLinked<RepcaExternalVariables::IBRANCHR>(),
            "ibranchr");
        check_required_signal(
            signals_.template isAttached<RepcaExternalVariables::IBRANCHI>(),
            signals_.template isAttached<RepcaExternalVariables::IBRANCHI>()
                && signals_.template isLinked<RepcaExternalVariables::IBRANCHI>(),
            "ibranchi");
        check_required_signal(
            signals_.template isAttached<RepcaExternalVariables::PBRANCH>(),
            signals_.template isAttached<RepcaExternalVariables::PBRANCH>()
                && signals_.template isLinked<RepcaExternalVariables::PBRANCH>(),
            "pbranch");
        check_required_signal(
            signals_.template isAttached<RepcaExternalVariables::QBRANCH>(),
            signals_.template isAttached<RepcaExternalVariables::QBRANCH>()
                && signals_.template isLinked<RepcaExternalVariables::QBRANCH>(),
            "qbranch");
        check_required_signal(
            signals_.template isAttached<RepcaExternalVariables::FREQ>(),
            signals_.template isAttached<RepcaExternalVariables::FREQ>()
                && signals_.template isLinked<RepcaExternalVariables::FREQ>(),
            "freq");

        auto check_optional_signal = [&](bool attached, bool linked, const char* name)
        {
          if (attached && !linked)
          {
            Log::error() << "Repca: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_optional_signal(
            signals_.template isAttached<RepcaExternalVariables::FREQREF>(),
            signals_.template isAttached<RepcaExternalVariables::FREQREF>()
                && signals_.template isLinked<RepcaExternalVariables::FREQREF>(),
            "freqref");
        check_optional_signal(
            signals_.template isAttached<RepcaExternalVariables::VREF>(),
            signals_.template isAttached<RepcaExternalVariables::VREF>()
                && signals_.template isLinked<RepcaExternalVariables::VREF>(),
            "vref");
        check_optional_signal(
            signals_.template isAttached<RepcaExternalVariables::QREF>(),
            signals_.template isAttached<RepcaExternalVariables::QREF>()
                && signals_.template isLinked<RepcaExternalVariables::QREF>(),
            "qref");
        check_optional_signal(
            signals_.template isAttached<RepcaExternalVariables::PPLANTREF>(),
            signals_.template isAttached<RepcaExternalVariables::PPLANTREF>()
                && signals_.template isLinked<RepcaExternalVariables::PPLANTREF>(),
            "pplantref");

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::initialize()
      {
        if (verify() > 0)
        {
          Log::error() << "Repca: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

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

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();
        const ScalarT ibranchr =
            signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHR>();
        const ScalarT ibranchi =
            signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHI>();
        const ScalarT pbranch =
            signals_.template readExternalVariable<RepcaExternalVariables::PBRANCH>();
        const ScalarT qbranch =
            signals_.template readExternalVariable<RepcaExternalVariables::QBRANCH>();
        const ScalarT freq =
            signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();

        const ScalarT qext0 = toComponentBase(y[QEXT]);
        const ScalarT pext0 = toComponentBase(y[PEXT]);

        const ScalarT vldc_r = vr - Rc_ * ibranchr + Xc_ * ibranchi;
        const ScalarT vldc_i = vi - Rc_ * ibranchi - Xc_ * ibranchr;

        y[V]      = std::sqrt(vr * vr + vi * vi);
        y[VLDC]   = std::sqrt(vldc_r * vldc_r + vldc_i * vldc_i);
        y[VDROOP] = y[V] + Kc_ * toComponentBase(qbranch);
        y[VCTRL]  = vcomp_on_ * y[VLDC] + vcomp_off_ * y[VDROOP];

        y[EF]                = ZERO<RealT>;
        const ScalarT pfreq0 = Ddn_ * Math::ramp(y[EF]) - Dup_ * Math::ramp(-y[EF]);

        y[VMEAS] = y[VCTRL];
        y[QMEAS] = toComponentBase(qbranch);
        y[PMEAS] = toComponentBase(pbranch);
        y[SFRZ]  = Math::above(y[V], Vfrz_);

        y[ERQ]    = ZERO<RealT>;
        y[ERQDB]  = Math::deadband2(y[ERQ], dbdlow_, dbdupper_);
        y[ERQLIM] = Math::clamp(y[ERQDB], emin_, emax_);
        y[QPI]    = qext0;
        y[XQLAG]  = qext0;
        y[QEXT]   = toSystemBase(qext0);
        y[XQPI]   = qext0 - Kp_ * y[ERQLIM];

        y[EP]    = ZERO<RealT>;
        y[EPLIM] = Math::clamp(y[EP], femin_, femax_);
        y[PREF]  = Freqflag_ ? pext0 : Math::clamp(y[PMEAS], Pmin_, Pmax_);
        y[PPI]   = y[PREF];
        y[XPPI]  = y[PREF] - Kpg_ * y[EPLIM];
        y[PEXT]  = toSystemBase(freq_on_ * y[PREF]);

        const ScalarT q_aw = Math::antiwindup(y[QPI], Ki_ * y[ERQLIM], Qmin_, Qmax_);
        const ScalarT p_aw = Math::antiwindup(y[PPI], Kig_ * y[EPLIM], Pmin_, Pmax_);
        if (std::abs(static_cast<RealT>(q_aw)) > static_cast<RealT>(1.0e-10))
        {
          Log::error() << "Repca: reactive PI antiwindup rate is nonzero at initialization\n";
          return 1;
        }
        if (std::abs(static_cast<RealT>(p_aw)) > static_cast<RealT>(1.0e-10))
        {
          Log::error() << "Repca: active PI antiwindup rate is nonzero at initialization\n";
          return 1;
        }

        freqref_set_   = freq;
        vref_set_      = y[VMEAS];
        qref_set_      = toSystemBase(y[QMEAS]);
        pplantref_set_ = toSystemBase(y[PMEAS] - pfreq0);

        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::FREQREF>(
              freqref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::VREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::VREF>(vref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::QREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::QREF>(qref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::PPLANTREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::PPLANTREF>(
              pplantref_set_);
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
      int Repca<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(RepcaInternalVariables::VMEAS)] = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::QMEAS)] = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::XQPI)]  = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::XQLAG)] = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::PMEAS)] = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::XPPI)]  = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::PREF)]  = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int Repca<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
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

        const auto IBRANCHR  = static_cast<size_t>(RepcaExternalVariables::IBRANCHR);
        const auto IBRANCHI  = static_cast<size_t>(RepcaExternalVariables::IBRANCHI);
        const auto PBRANCH   = static_cast<size_t>(RepcaExternalVariables::PBRANCH);
        const auto QBRANCH   = static_cast<size_t>(RepcaExternalVariables::QBRANCH);
        const auto FREQ      = static_cast<size_t>(RepcaExternalVariables::FREQ);
        const auto FREQREF   = static_cast<size_t>(RepcaExternalVariables::FREQREF);
        const auto VREF      = static_cast<size_t>(RepcaExternalVariables::VREF);
        const auto QREF      = static_cast<size_t>(RepcaExternalVariables::QREF);
        const auto PPLANTREF = static_cast<size_t>(RepcaExternalVariables::PPLANTREF);

        const ScalarT vmeas  = y[VMEAS];
        const ScalarT qmeas  = y[QMEAS];
        const ScalarT xqpi   = y[XQPI];
        const ScalarT xqlag  = y[XQLAG];
        const ScalarT pmeas  = y[PMEAS];
        const ScalarT xppi   = y[XPPI];
        const ScalarT pref   = y[PREF];
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
        const ScalarT pref_dot  = yp[PREF];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT ibranchr  = ws[IBRANCHR];
        const ScalarT ibranchi  = ws[IBRANCHI];
        const ScalarT pbranch   = ws[PBRANCH];
        const ScalarT qbranch   = ws[QBRANCH];
        const ScalarT freq      = ws[FREQ];
        const ScalarT freqref   = ws[FREQREF];
        const ScalarT vref      = ws[VREF];
        const ScalarT qref      = toComponentBase(ws[QREF]);
        const ScalarT pplantref = toComponentBase(ws[PPLANTREF]);

        const ScalarT vldc_r = vr - Rc_ * ibranchr + Xc_ * ibranchi;
        const ScalarT vldc_i = vi - Rc_ * ibranchi - Xc_ * ibranchr;
        const ScalarT pfreq  = Ddn_ * Math::ramp(ef) - Dup_ * Math::ramp(-ef);

        f[VMEAS] = -vmeas_dot + (vctrl - vmeas) / Tfltr_;
        f[QMEAS] = -qmeas_dot + (toComponentBase(qbranch) - qmeas) / Tfltr_;
        f[XQPI]  = -xqpi_dot + sfrz * Math::antiwindup(qpi, Ki_ * erqlim, Qmin_, Qmax_);
        f[XQLAG] = -xqlag_dot + (qpi - xqlag) / Tfv_;
        f[PMEAS] = -pmeas_dot + (toComponentBase(pbranch) - pmeas) / Tp_;
        f[XPPI]  = -xppi_dot + Math::antiwindup(ppi, Kig_ * eplim, Pmin_, Pmax_);
        f[PREF]  = -pref_dot + (ppi - pref) / Tlag_;

        f[V]      = -v * v + vr * vr + vi * vi;
        f[VLDC]   = -vldc * vldc + vldc_r * vldc_r + vldc_i * vldc_i;
        f[VDROOP] = -vdroop + v + Kc_ * toComponentBase(qbranch);
        f[VCTRL]  = -vctrl + vcomp_on_ * vldc + vcomp_off_ * vdroop;
        f[SFRZ]   = -sfrz + Math::above(v, Vfrz_);
        f[ERQ]    = -erq + ref_on_ * (vref - vmeas) + ref_off_ * (qref - qmeas);
        f[ERQDB]  = -erqdb + Math::deadband2(erq, dbdlow_, dbdupper_);
        f[ERQLIM] = -erqlim + Math::clamp(erqdb, emin_, emax_);
        f[QPI]    = -qpi + Math::clamp(Kp_ * erqlim + xqpi, Qmin_, Qmax_);
        f[QEXT]   = -Tfv_ * (qext - xqlag) + Tft_ * (qpi - xqlag);

        f[EF]    = -ef + Math::deadband2(freqref - freq, fdbd1_, fdbd2_);
        f[EP]    = -ep + pplantref - pmeas + pfreq;
        f[EPLIM] = -eplim + Math::clamp(ep, femin_, femax_);
        f[PPI]   = -ppi + Math::clamp(Kpg_ * eplim + xppi, Pmin_, Pmax_);
        f[PEXT]  = -pext + freq_on_ * pref;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateResidual()
      {
        const auto IBRANCHR  = static_cast<size_t>(RepcaExternalVariables::IBRANCHR);
        const auto IBRANCHI  = static_cast<size_t>(RepcaExternalVariables::IBRANCHI);
        const auto PBRANCH   = static_cast<size_t>(RepcaExternalVariables::PBRANCH);
        const auto QBRANCH   = static_cast<size_t>(RepcaExternalVariables::QBRANCH);
        const auto FREQ      = static_cast<size_t>(RepcaExternalVariables::FREQ);
        const auto FREQREF   = static_cast<size_t>(RepcaExternalVariables::FREQREF);
        const auto VREF      = static_cast<size_t>(RepcaExternalVariables::VREF);
        const auto QREF      = static_cast<size_t>(RepcaExternalVariables::QREF);
        const auto PPLANTREF = static_cast<size_t>(RepcaExternalVariables::PPLANTREF);

        std::fill(ws_.begin(), ws_.end(), ZERO<RealT>);
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        ws_[FREQREF]   = freqref_set_;
        ws_[VREF]      = vref_set_;
        ws_[QREF]      = qref_set_;
        ws_[PPLANTREF] = pplantref_set_;

        if (signals_.template isAttached<RepcaExternalVariables::IBRANCHR>())
        {
          ws_[IBRANCHR] =
              signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHR>();
          ws_indices_[IBRANCHR] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::IBRANCHR>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::IBRANCHI>())
        {
          ws_[IBRANCHI] =
              signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHI>();
          ws_indices_[IBRANCHI] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::IBRANCHI>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::PBRANCH>())
        {
          ws_[PBRANCH] =
              signals_.template readExternalVariable<RepcaExternalVariables::PBRANCH>();
          ws_indices_[PBRANCH] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::PBRANCH>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::QBRANCH>())
        {
          ws_[QBRANCH] =
              signals_.template readExternalVariable<RepcaExternalVariables::QBRANCH>();
          ws_indices_[QBRANCH] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::QBRANCH>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::FREQ>())
        {
          ws_[FREQ] = signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();
          ws_indices_[FREQ] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQ>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
        {
          ws_[FREQREF] =
              signals_.template readExternalVariable<RepcaExternalVariables::FREQREF>();
          ws_indices_[FREQREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::VREF>())
        {
          ws_[VREF] = signals_.template readExternalVariable<RepcaExternalVariables::VREF>();
          ws_indices_[VREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::VREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::QREF>())
        {
          ws_[QREF] = signals_.template readExternalVariable<RepcaExternalVariables::QREF>();
          ws_indices_[QREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::QREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::PPLANTREF>())
        {
          ws_[PPLANTREF] =
              signals_.template readExternalVariable<RepcaExternalVariables::PPLANTREF>();
          ws_indices_[PPLANTREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::PPLANTREF>();
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

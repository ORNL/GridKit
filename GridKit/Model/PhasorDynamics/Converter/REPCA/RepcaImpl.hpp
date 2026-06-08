/**
 * @file RepcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REPCA plant-control model.
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
      }

      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(RepcaInternalVariables::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::~Repca()
      {
      }

      template <typename scalar_type, typename index_type>
      scalar_type& Repca<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      template <typename scalar_type, typename index_type>
      scalar_type& Repca<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::setDerivedParameters()
      {
        system_to_component_base_ = ONE<RealT>;
        if (mva_base_ > ZERO<RealT>)
        {
          const RealT va_component_base = mva_base_ * static_cast<RealT>(1.0e6);
          system_to_component_base_     = va_system_base_ / va_component_base;
        }
        vcomp_off_ = ONE<RealT> - VcompFlag_;
        ref_off_   = ONE<RealT> - RefFlag_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Repca<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * system_to_component_base_;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::verifyBranchSignalPorts() const
      {
        int ret = 0;

        const bool has_pbranch  = signals_.template isAttached<RepcaExternalVariables::PBRANCH>();
        const bool has_qbranch  = signals_.template isAttached<RepcaExternalVariables::QBRANCH>();
        const bool has_ibranchr = signals_.template isAttached<RepcaExternalVariables::IBRANCHR>();
        const bool has_ibranchi = signals_.template isAttached<RepcaExternalVariables::IBRANCHI>();

        if (!has_pbranch || !has_qbranch)
        {
          Log::error() << "Repca: pbranch and qbranch must be connected\n";
          ret += 1;
        }

        if (!has_ibranchr || !has_ibranchi)
        {
          Log::error() << "Repca: ibranchr and ibranchi must be connected\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        auto load_required_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Repca: missing required parameter '" << name << "'\n";
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
            Log::error() << "Repca: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto load_required_switch = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Repca: missing required parameter '" << name << "'\n";
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
            Log::error() << "Repca: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        load_required_real(Params::mva, mva_base_, "mva");
        load_required_switch(Params::VcompFlag, VcompFlag_, "VcompFlag");
        load_required_switch(Params::RefFlag, RefFlag_, "RefFlag");
        load_required_switch(Params::Freqflag, Freqflag_, "Freqflag");
        load_required_real(Params::Tfltr, Tfltr_, "Tfltr");
        load_required_real(Params::Tft, Tft_, "Tft");
        load_required_real(Params::Tfv, Tfv_, "Tfv");
        load_required_real(Params::Tp, Tp_, "Tp");
        load_required_real(Params::Tlag, Tlag_, "Tlag");
        load_required_real(Params::Vfrz, Vfrz_, "Vfrz");
        load_required_real(Params::Rc, Rc_, "Rc");
        load_required_real(Params::Xc, Xc_, "Xc");
        load_required_real(Params::Kc, Kc_, "Kc");
        load_required_real(Params::dbdlow, dbdlow_, "dbdlow");
        load_required_real(Params::dbdupper, dbdupper_, "dbdupper");
        load_required_real(Params::emax, emax_, "emax");
        load_required_real(Params::emin, emin_, "emin");
        load_required_real(Params::Kp, Kp_, "Kp");
        load_required_real(Params::Ki, Ki_, "Ki");
        load_required_real(Params::Qmax, Qmax_, "Qmax");
        load_required_real(Params::Qmin, Qmin_, "Qmin");
        load_required_real(Params::fdbd1, fdbd1_, "fdbd1");
        load_required_real(Params::fdbd2, fdbd2_, "fdbd2");
        load_required_real(Params::Ddn, Ddn_, "Ddn");
        load_required_real(Params::Dup, Dup_, "Dup");
        load_required_real(Params::femax, femax_, "femax");
        load_required_real(Params::femin, femin_, "femin");
        load_required_real(Params::Kpg, Kpg_, "Kpg");
        load_required_real(Params::Kig, Kig_, "Kig");
        load_required_real(Params::Pmax, Pmax_, "Pmax");
        load_required_real(Params::Pmin, Pmin_, "Pmin");

        setDerivedParameters();
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
                      { return y_[index(RepcaInternalVariables::QEXT)]; });
        monitor_->set(Variable::pext, [this, index]
                      { return y_[index(RepcaInternalVariables::PEXT)]; });
        monitor_->set(Variable::vmeas, [this, index]
                      { return y_[index(RepcaInternalVariables::VMEAS)]; });
        monitor_->set(Variable::qmeas, [this, index]
                      { return y_[index(RepcaInternalVariables::QMEAS)]; });
        monitor_->set(Variable::pmeas, [this, index]
                      { return y_[index(RepcaInternalVariables::PMEAS)]; });
        monitor_->set(Variable::pref, [this, index]
                      { return y_[index(RepcaInternalVariables::PREF)]; });
        monitor_->set(Variable::vctrl, [this, index]
                      { return y_[index(RepcaInternalVariables::VCTRL)]; });
        monitor_->set(Variable::sfrz, [this, index]
                      { return y_[index(RepcaInternalVariables::SFRZ)]; });
        monitor_->set(Variable::qpi, [this, index]
                      { return y_[index(RepcaInternalVariables::QPI)]; });
        monitor_->set(Variable::pfreq, [this, index]
                      {
                        const auto EF = index(RepcaInternalVariables::EF);
                        return Ddn_ * Math::ramp(y_[EF]) - Dup_ * Math::ramp(-y_[EF]); });
        monitor_->set(Variable::ppi, [this, index]
                      { return y_[index(RepcaInternalVariables::PPI)]; });
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
        setDerivedParameters();

        size_     = static_cast<IdxT>(RepcaInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        f_.assign(size, ScalarT{0});
        y_.assign(size, ScalarT{0});
        yp_.assign(size, ScalarT{0});
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
          signals_.template getSignalNode<RepcaInternalVariables::QEXT>()->set(
              &y_[static_cast<size_t>(RepcaInternalVariables::QEXT)],
              &(this->getVariableIndex(static_cast<IdxT>(RepcaInternalVariables::QEXT))));
        }

        if (signals_.template isAssigned<RepcaInternalVariables::PEXT>())
        {
          signals_.template getSignalNode<RepcaInternalVariables::PEXT>()->set(
              &y_[static_cast<size_t>(RepcaInternalVariables::PEXT)],
              &(this->getVariableIndex(static_cast<IdxT>(RepcaInternalVariables::PEXT))));
        }

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

        check(mva_base_ >= ZERO<RealT>, "mva must be non-negative");
        check(VcompFlag_ == ZERO<RealT> || VcompFlag_ == ONE<RealT>,
              "VcompFlag must be 0 or 1");
        check(RefFlag_ == ZERO<RealT> || RefFlag_ == ONE<RealT>, "RefFlag must be 0 or 1");
        check(Freqflag_ == ZERO<RealT> || Freqflag_ == ONE<RealT>,
              "Freqflag must be 0 or 1");
        check(Tfltr_ >= ZERO<RealT>, "Tfltr must be non-negative");
        check(Tp_ >= ZERO<RealT>, "Tp must be non-negative");
        check(Tlag_ >= ZERO<RealT>, "Tlag must be non-negative");
        check(Tft_ >= ZERO<RealT>, "Tft must be non-negative");
        check(Tfv_ > ZERO<RealT>, "Tfv must be positive");
        check(Vfrz_ >= ZERO<RealT>, "Vfrz must be non-negative");
        check(dbdlow_ <= ZERO<RealT> && ZERO<RealT> <= dbdupper_,
              "dbdlow <= 0 <= dbdupper is required");
        check(emin_ <= ZERO<RealT> && ZERO<RealT> <= emax_, "emin <= 0 <= emax is required");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(fdbd1_ <= ZERO<RealT> && ZERO<RealT> <= fdbd2_,
              "fdbd1 <= 0 <= fdbd2 is required");
        check(Ddn_ >= ZERO<RealT>, "Ddn must be non-negative");
        check(Dup_ >= ZERO<RealT>, "Dup must be non-negative");
        check(femin_ <= ZERO<RealT> && ZERO<RealT> <= femax_, "femin <= 0 <= femax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");
        ret += verifyBranchSignalPorts();

        if (signals_.template isAttached<RepcaExternalVariables::IBRANCHR>()
            && !signals_.template isLinked<RepcaExternalVariables::IBRANCHR>())
        {
          Log::error() << "Repca: ibranchr signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAttached<RepcaExternalVariables::IBRANCHI>()
            && !signals_.template isLinked<RepcaExternalVariables::IBRANCHI>())
        {
          Log::error() << "Repca: ibranchi signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAttached<RepcaExternalVariables::QBRANCH>()
            && !signals_.template isLinked<RepcaExternalVariables::QBRANCH>())
        {
          Log::error() << "Repca: qbranch signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAttached<RepcaExternalVariables::PBRANCH>()
            && !signals_.template isLinked<RepcaExternalVariables::PBRANCH>())
        {
          Log::error() << "Repca: pbranch signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAttached<RepcaExternalVariables::VREF>()
            && !signals_.template isLinked<RepcaExternalVariables::VREF>())
        {
          Log::error() << "Repca: vref signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAttached<RepcaExternalVariables::QREF>()
            && !signals_.template isLinked<RepcaExternalVariables::QREF>())
        {
          Log::error() << "Repca: qref signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAttached<RepcaExternalVariables::PPLANTREF>()
            && !signals_.template isLinked<RepcaExternalVariables::PPLANTREF>())
        {
          Log::error() << "Repca: pplantref signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAttached<RepcaExternalVariables::FREQ>()
            && !signals_.template isLinked<RepcaExternalVariables::FREQ>())
        {
          Log::error() << "Repca: freq signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>()
            && !signals_.template isLinked<RepcaExternalVariables::FREQREF>())
        {
          Log::error() << "Repca: freqref signal attached with no linked source\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Repca<scalar_type, index_type>::readExternalOrDefault(
          RepcaExternalVariables variable,
          scalar_type            default_value) const
      {
        switch (variable)
        {
        case RepcaExternalVariables::IBRANCHR:
          if (signals_.template isAttached<RepcaExternalVariables::IBRANCHR>())
            return signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHR>();
          break;
        case RepcaExternalVariables::IBRANCHI:
          if (signals_.template isAttached<RepcaExternalVariables::IBRANCHI>())
            return signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHI>();
          break;
        case RepcaExternalVariables::QBRANCH:
          if (signals_.template isAttached<RepcaExternalVariables::QBRANCH>())
            return signals_.template readExternalVariable<RepcaExternalVariables::QBRANCH>();
          break;
        case RepcaExternalVariables::PBRANCH:
          if (signals_.template isAttached<RepcaExternalVariables::PBRANCH>())
            return signals_.template readExternalVariable<RepcaExternalVariables::PBRANCH>();
          break;
        case RepcaExternalVariables::VREF:
          if (signals_.template isAttached<RepcaExternalVariables::VREF>())
            return signals_.template readExternalVariable<RepcaExternalVariables::VREF>();
          break;
        case RepcaExternalVariables::QREF:
          if (signals_.template isAttached<RepcaExternalVariables::QREF>())
            return signals_.template readExternalVariable<RepcaExternalVariables::QREF>();
          break;
        case RepcaExternalVariables::PPLANTREF:
          if (signals_.template isAttached<RepcaExternalVariables::PPLANTREF>())
            return signals_.template readExternalVariable<RepcaExternalVariables::PPLANTREF>();
          break;
        case RepcaExternalVariables::FREQ:
          if (signals_.template isAttached<RepcaExternalVariables::FREQ>())
            return signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();
          break;
        case RepcaExternalVariables::FREQREF:
          if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
            return signals_.template readExternalVariable<RepcaExternalVariables::FREQREF>();
          break;
        case RepcaExternalVariables::MAXIMUM:
          break;
        }
        return default_value;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::initialize()
      {
        if (parameter_error_count_ > 0 || verify() > 0)
        {
          Log::error() << "Repca: cannot initialize with invalid configuration\n";
          return 1;
        }

        setDerivedParameters();

        const auto VMEAS  = static_cast<size_t>(RepcaInternalVariables::VMEAS);
        const auto QMEAS  = static_cast<size_t>(RepcaInternalVariables::QMEAS);
        const auto XQ     = static_cast<size_t>(RepcaInternalVariables::XQ);
        const auto XQEXT  = static_cast<size_t>(RepcaInternalVariables::XQEXT);
        const auto PMEAS  = static_cast<size_t>(RepcaInternalVariables::PMEAS);
        const auto XP     = static_cast<size_t>(RepcaInternalVariables::XP);
        const auto PREF   = static_cast<size_t>(RepcaInternalVariables::PREF);
        const auto VREG   = static_cast<size_t>(RepcaInternalVariables::VREG);
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

        const ScalarT vr    = Vr();
        const ScalarT vi    = Vi();
        const ScalarT vt_sq = vr * vr + vi * vi;

        if (static_cast<RealT>(vt_sq) <= ZERO<RealT>)
        {
          Log::error()
              << "Repca: regulated-bus voltage magnitude must be positive at initialization\n";
          return 1;
        }

        const ScalarT pbr0_system =
            signals_.template readExternalVariable<RepcaExternalVariables::PBRANCH>();
        const ScalarT qbr0_system =
            signals_.template readExternalVariable<RepcaExternalVariables::QBRANCH>();
        const ScalarT pbr0 = toComponentBase(pbr0_system);
        const ScalarT qbr0 = toComponentBase(qbr0_system);
        const ScalarT ibrr0 =
            signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHR>();
        const ScalarT ibri0 =
            signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHI>();

        y_[VREG]   = std::sqrt(vr * vr + vi * vi);
        y_[VLDC]   = std::sqrt((vr - Rc_ * ibrr0 + Xc_ * ibri0)
                                 * (vr - Rc_ * ibrr0 + Xc_ * ibri0)
                             + (vi - Rc_ * ibri0 - Xc_ * ibrr0)
                                   * (vi - Rc_ * ibri0 - Xc_ * ibrr0));
        y_[VDROOP] = y_[VREG] + Kc_ * qbr0;
        y_[VCTRL]  = VcompFlag_ * y_[VLDC] + vcomp_off_ * y_[VDROOP];

        y_[VMEAS] = y_[VCTRL];
        y_[QMEAS] = qbr0;
        y_[PMEAS] = pbr0;

        freq_set_    = readExternalOrDefault(RepcaExternalVariables::FREQ, ScalarT{ZERO<RealT>});
        freqref_set_ = readExternalOrDefault(RepcaExternalVariables::FREQREF, ScalarT{ZERO<RealT>});

        y_[EF]               = Math::deadband2(freqref_set_ - freq_set_, fdbd1_, fdbd2_);
        const ScalarT pfreq0 = Ddn_ * Math::ramp(y_[EF]) - Dup_ * Math::ramp(-y_[EF]);

        vref_set_ = readExternalOrDefault(RepcaExternalVariables::VREF, y_[VMEAS]);
        qref_set_ = readExternalOrDefault(RepcaExternalVariables::QREF, y_[QMEAS]);
        pplantref_set_ =
            readExternalOrDefault(RepcaExternalVariables::PPLANTREF, y_[PMEAS] - pfreq0);

        y_[SFRZ]   = Math::above(y_[VREG], Vfrz_);
        y_[ERQ]    = RefFlag_ * (vref_set_ - y_[VMEAS]) + ref_off_ * (qref_set_ - y_[QMEAS]);
        y_[ERQDB]  = Math::deadband2(y_[ERQ], dbdlow_, dbdupper_);
        y_[ERQLIM] = Math::clamp(y_[ERQDB], emin_, emax_);

        ScalarT qext0 = Math::clamp(qbr0, Qmin_, Qmax_);
        if (signals_.template isAssigned<RepcaInternalVariables::QEXT>())
        {
          qext0 = y_[QEXT];
        }

        y_[QEXT]  = qext0;
        y_[QPI]   = y_[QEXT];
        y_[XQEXT] = y_[QEXT];
        y_[XQ]    = y_[QEXT] - Kp_ * y_[ERQLIM];

        y_[EP]    = pplantref_set_ - y_[PMEAS] + pfreq0;
        y_[EPLIM] = Math::clamp(y_[EP], femin_, femax_);

        ScalarT pref0 = Math::clamp(pbr0, Pmin_, Pmax_);
        if (signals_.template isAssigned<RepcaInternalVariables::PEXT>()
            && Freqflag_ != ZERO<RealT>)
        {
          pref0 = y_[PEXT] / Freqflag_;
        }

        y_[PREF] = pref0;
        y_[PPI]  = y_[PREF];
        y_[XP]   = y_[PREF] - Kpg_ * y_[EPLIM];
        y_[PEXT] = Freqflag_ * y_[PREF];

        static constexpr RealT init_tol = static_cast<RealT>(1.0e-10);

        std::fill(yp_.begin(), yp_.end(), ZERO<RealT>);

        const ScalarT q_rate =
            y_[SFRZ] * Math::antiwindup(y_[QPI], Ki_ * y_[ERQLIM], Qmin_, Qmax_);

        const ScalarT p_rate =
            Math::antiwindup(y_[PPI], Kig_ * y_[EPLIM], Pmin_, Pmax_);

        if (std::abs(static_cast<RealT>(q_rate)) > init_tol)
        {
          Log::error() << "Repca: reactive controller initial condition is not steady state\n";
          return 1;
        }

        if (std::abs(static_cast<RealT>(p_rate)) > init_tol)
        {
          Log::error() << "Repca: active-power controller initial condition is not steady state\n";
          return 1;
        }

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(RepcaInternalVariables::VMEAS)] = (Tfltr_ > ZERO<RealT>);
        tag_[static_cast<size_t>(RepcaInternalVariables::QMEAS)] = (Tfltr_ > ZERO<RealT>);
        tag_[static_cast<size_t>(RepcaInternalVariables::XQ)]    = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::XQEXT)] = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::PMEAS)] = (Tp_ > ZERO<RealT>);
        tag_[static_cast<size_t>(RepcaInternalVariables::XP)]    = true;
        tag_[static_cast<size_t>(RepcaInternalVariables::PREF)]  = (Tlag_ > ZERO<RealT>);
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Repca<scalar_type, index_type>::evaluateInternalResidual(
          ScalarT* y,
          ScalarT* yp,
          ScalarT* wb,
          ScalarT* ws,
          ScalarT* f)
      {
        const auto VMEAS  = static_cast<size_t>(RepcaInternalVariables::VMEAS);
        const auto QMEAS  = static_cast<size_t>(RepcaInternalVariables::QMEAS);
        const auto XQ     = static_cast<size_t>(RepcaInternalVariables::XQ);
        const auto XQEXT  = static_cast<size_t>(RepcaInternalVariables::XQEXT);
        const auto PMEAS  = static_cast<size_t>(RepcaInternalVariables::PMEAS);
        const auto XP     = static_cast<size_t>(RepcaInternalVariables::XP);
        const auto PREF   = static_cast<size_t>(RepcaInternalVariables::PREF);
        const auto VREG   = static_cast<size_t>(RepcaInternalVariables::VREG);
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
        const auto QBRANCH   = static_cast<size_t>(RepcaExternalVariables::QBRANCH);
        const auto PBRANCH   = static_cast<size_t>(RepcaExternalVariables::PBRANCH);
        const auto VREF      = static_cast<size_t>(RepcaExternalVariables::VREF);
        const auto QREF      = static_cast<size_t>(RepcaExternalVariables::QREF);
        const auto PPLANTREF = static_cast<size_t>(RepcaExternalVariables::PPLANTREF);
        const auto FREQ      = static_cast<size_t>(RepcaExternalVariables::FREQ);
        const auto FREQREF   = static_cast<size_t>(RepcaExternalVariables::FREQREF);

        const ScalarT vr   = wb[0];
        const ScalarT vi   = wb[1];
        const ScalarT pbr  = toComponentBase(ws[PBRANCH]);
        const ScalarT qbr  = toComponentBase(ws[QBRANCH]);
        const ScalarT ibrr = ws[IBRANCHR];
        const ScalarT ibri = ws[IBRANCHI];

        f[VMEAS] = -Tfltr_ * yp[VMEAS] - y[VMEAS] + y[VCTRL];
        f[QMEAS] = -Tfltr_ * yp[QMEAS] - y[QMEAS] + qbr;
        f[XQ]    = -yp[XQ]
                + y[SFRZ] * Math::antiwindup(y[QPI], Ki_ * y[ERQLIM], Qmin_, Qmax_);
        f[XQEXT] = -Tfv_ * yp[XQEXT] - y[XQEXT] + y[QPI];
        f[PMEAS] = -Tp_ * yp[PMEAS] - y[PMEAS] + pbr;
        f[XP]    = -yp[XP] + Math::antiwindup(y[PPI], Kig_ * y[EPLIM], Pmin_, Pmax_);
        f[PREF]  = -Tlag_ * yp[PREF] - y[PREF] + y[PPI];

        f[VREG] = -y[VREG] * y[VREG] + vr * vr + vi * vi;
        f[VLDC] = -y[VLDC] * y[VLDC]
                  + (vr - Rc_ * ibrr + Xc_ * ibri)
                        * (vr - Rc_ * ibrr + Xc_ * ibri)
                  + (vi - Rc_ * ibri - Xc_ * ibrr)
                        * (vi - Rc_ * ibri - Xc_ * ibrr);
        f[VDROOP] = -y[VDROOP] + y[VREG] + Kc_ * qbr;
        f[VCTRL]  = -y[VCTRL] + VcompFlag_ * y[VLDC] + vcomp_off_ * y[VDROOP];
        f[SFRZ]   = -y[SFRZ] + Math::above(y[VREG], Vfrz_);
        f[ERQ]    = -y[ERQ] + RefFlag_ * (ws[VREF] - y[VMEAS])
                 + ref_off_ * (ws[QREF] - y[QMEAS]);
        f[ERQDB]            = -y[ERQDB] + Math::deadband2(y[ERQ], dbdlow_, dbdupper_);
        f[ERQLIM]           = -y[ERQLIM] + Math::clamp(y[ERQDB], emin_, emax_);
        f[QPI]              = -y[QPI] + Math::clamp(Kp_ * y[ERQLIM] + y[XQ], Qmin_, Qmax_);
        f[QEXT]             = -Tfv_ * (y[QEXT] - y[XQEXT]) + Tft_ * (y[QPI] - y[XQEXT]);
        f[EF]               = -y[EF] + Math::deadband2(ws[FREQREF] - ws[FREQ], fdbd1_, fdbd2_);
        const ScalarT pfreq = Ddn_ * Math::ramp(y[EF]) - Dup_ * Math::ramp(-y[EF]);
        f[EP]               = -y[EP] + ws[PPLANTREF] - y[PMEAS] + pfreq;
        f[EPLIM]            = -y[EPLIM] + Math::clamp(y[EP], femin_, femax_);
        f[PPI]              = -y[PPI] + Math::clamp(Kpg_ * y[EPLIM] + y[XP], Pmin_, Pmax_);
        f[PEXT]             = -y[PEXT] + Freqflag_ * y[PREF];

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateResidual()
      {
        setDerivedParameters();

        const auto IBRANCHR  = static_cast<size_t>(RepcaExternalVariables::IBRANCHR);
        const auto IBRANCHI  = static_cast<size_t>(RepcaExternalVariables::IBRANCHI);
        const auto QBRANCH   = static_cast<size_t>(RepcaExternalVariables::QBRANCH);
        const auto PBRANCH   = static_cast<size_t>(RepcaExternalVariables::PBRANCH);
        const auto VREF      = static_cast<size_t>(RepcaExternalVariables::VREF);
        const auto QREF      = static_cast<size_t>(RepcaExternalVariables::QREF);
        const auto PPLANTREF = static_cast<size_t>(RepcaExternalVariables::PPLANTREF);
        const auto FREQ      = static_cast<size_t>(RepcaExternalVariables::FREQ);
        const auto FREQREF   = static_cast<size_t>(RepcaExternalVariables::FREQREF);

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();

        ws_[IBRANCHR]  = ZERO<RealT>;
        ws_[IBRANCHI]  = ZERO<RealT>;
        ws_[QBRANCH]   = ZERO<RealT>;
        ws_[PBRANCH]   = ZERO<RealT>;
        ws_[VREF]      = vref_set_;
        ws_[QREF]      = qref_set_;
        ws_[PPLANTREF] = pplantref_set_;
        ws_[FREQ]      = freq_set_;
        ws_[FREQREF]   = freqref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

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
        if (signals_.template isAttached<RepcaExternalVariables::QBRANCH>())
        {
          ws_[QBRANCH] = signals_.template readExternalVariable<RepcaExternalVariables::QBRANCH>();
          ws_indices_[QBRANCH] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::QBRANCH>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::PBRANCH>())
        {
          ws_[PBRANCH] = signals_.template readExternalVariable<RepcaExternalVariables::PBRANCH>();
          ws_indices_[PBRANCH] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::PBRANCH>();
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
        if (signals_.template isAttached<RepcaExternalVariables::FREQ>())
        {
          ws_[FREQ] = signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();
          ws_indices_[FREQ] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQ>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
        {
          ws_[FREQREF] = signals_.template readExternalVariable<RepcaExternalVariables::FREQREF>();
          ws_indices_[FREQREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQREF>();
        }

        wb_[0] = vr;
        wb_[1] = vi;

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
        return 0;
      }
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

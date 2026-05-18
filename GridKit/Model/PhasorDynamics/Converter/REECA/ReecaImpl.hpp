/**
 * @file ReecaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REECA phasor-dynamics converter-control model.
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

      template <class ScalarT, typename IdxT>
      Reeca<ScalarT, IdxT>::Reeca(bus_type* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(ReecaInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Reeca<ScalarT, IdxT>::Reeca(bus_type* bus, const model_data_type& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(ReecaInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Reeca<ScalarT, IdxT>::~Reeca()
      {
      }

      template <class ScalarT, typename IdxT>
      void Reeca<ScalarT, IdxT>::setDerivedParameters()
      {
        spf_on_  = spf_ ? ONE<RealT> : ZERO<RealT>;
        spf_off_ = ONE<RealT> - spf_on_;
        sV_on_   = sV_ ? ONE<RealT> : ZERO<RealT>;
        sV_off_  = ONE<RealT> - sV_on_;
        sQ_on_   = sQ_ ? ONE<RealT> : ZERO<RealT>;
        sQ_off_  = ONE<RealT> - sQ_on_;
        sP_on_   = sP_ ? ONE<RealT> : ZERO<RealT>;
        sPQ_on_  = sPQ_ ? ONE<RealT> : ZERO<RealT>;
        sPQ_off_ = ONE<RealT> - sPQ_on_;
      }

      template <class ScalarT, typename IdxT>
      void Reeca<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
      {
        using External = ReecaExternalVariables;
        using Params   = typename model_data_type::Parameters;
        using Ports    = typename model_data_type::Ports;

        parameter_error_count_ = static_cast<IdxT>(data.input_error_count);
        has_Vref0_             = false;
        has_pe_init_           = false;
        has_qgen_init_         = false;
        has_omega_init_        = false;
        has_qext_init_         = false;
        has_pfaref_init_       = false;
        has_pref_init_         = false;
        pe_init_               = ZERO<RealT>;
        qgen_init_             = ZERO<RealT>;
        omega_init_            = ZERO<RealT>;
        qext_init_             = ZERO<RealT>;
        pfaref_init_           = ZERO<RealT>;
        pref_init_             = ZERO<RealT>;
        Thld2_                 = ZERO<RealT>;

        auto loadRequiredReal = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Reeca: missing required parameter '" << name << "'\n";
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
            Log::error() << "Reeca: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto loadOptionalReal = [&](auto key, RealT& target, const char* name, bool& present)
        {
          present = false;
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            target  = *real_value;
            present = true;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            target  = static_cast<RealT>(*index_value);
            present = true;
          }
          else
          {
            Log::error() << "Reeca: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto loadRequiredSwitch = [&](auto key, bool& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Reeca: missing required parameter '" << name << "'\n";
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
            Log::error() << "Reeca: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        for (const auto& [key, value] : data.external_initial_values)
        {
          static_cast<void>(value);
          if (key == External::MAXIMUM)
          {
            Log::error() << "Reeca: init key 'MAXIMUM' is not a valid external variable\n";
            ++parameter_error_count_;
          }
        }

        auto loadOptionalInitReal = [&](External key, RealT& target, const char* name, bool& present)
        {
          present = false;
          if (!data.external_initial_values.contains(key))
          {
            return;
          }

          const auto& value = data.external_initial_values.at(key);
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            target  = *real_value;
            present = true;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            target  = static_cast<RealT>(*index_value);
            present = true;
          }
          else
          {
            Log::error() << "Reeca: init value '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        loadRequiredReal(Params::Sbase, Sbase_, "Sbase");
        loadRequiredSwitch(Params::spf, spf_, "spf");
        loadRequiredSwitch(Params::sV, sV_, "sV");
        loadRequiredSwitch(Params::sQ, sQ_, "sQ");
        loadRequiredSwitch(Params::sP, sP_, "sP");
        loadRequiredSwitch(Params::sPQ, sPQ_, "sPQ");
        loadRequiredReal(Params::Trv, Trv_, "Trv");
        loadRequiredReal(Params::Tp, Tp_, "Tp");
        loadOptionalReal(Params::Vref0, Vref0_, "Vref0", has_Vref0_);
        loadRequiredReal(Params::Vdip, Vdip_, "Vdip");
        loadRequiredReal(Params::Vup, Vup_, "Vup");
        loadRequiredReal(Params::Dbd1, Dbd1_, "Dbd1");
        loadRequiredReal(Params::Dbd2, Dbd2_, "Dbd2");
        loadRequiredReal(Params::Kqv, Kqv_, "Kqv");
        loadRequiredReal(Params::Iql1, Iql1_, "Iql1");
        loadRequiredReal(Params::Iqh1, Iqh1_, "Iqh1");
        loadRequiredReal(Params::Iqfrz, Iqfrz_, "Iqfrz");
        loadRequiredReal(Params::Thld, Thld_, "Thld");
        loadRequiredReal(Params::Qmax, Qmax_, "Qmax");
        loadRequiredReal(Params::Qmin, Qmin_, "Qmin");
        loadRequiredReal(Params::Kqp, Kqp_, "Kqp");
        loadRequiredReal(Params::Kqi, Kqi_, "Kqi");
        loadRequiredReal(Params::Vmax, Vmax_, "Vmax");
        loadRequiredReal(Params::Vmin, Vmin_, "Vmin");
        loadRequiredReal(Params::Vref1, Vref1_, "Vref1");
        loadRequiredReal(Params::Kvp, Kvp_, "Kvp");
        loadRequiredReal(Params::Kvi, Kvi_, "Kvi");
        loadRequiredReal(Params::Tiq, Tiq_, "Tiq");
        loadRequiredReal(Params::Tpord, Tpord_, "Tpord");
        loadRequiredReal(Params::dPmax, dPmax_, "dPmax");
        loadRequiredReal(Params::dPmin, dPmin_, "dPmin");
        loadRequiredReal(Params::Pmax, Pmax_, "Pmax");
        loadRequiredReal(Params::Pmin, Pmin_, "Pmin");
        loadRequiredReal(Params::Imax, Imax_, "Imax");
        loadRequiredReal(Params::vq1, vq1_, "vq1");
        loadRequiredReal(Params::lq1, lq1_, "lq1");
        loadRequiredReal(Params::vq2, vq2_, "vq2");
        loadRequiredReal(Params::lq2, lq2_, "lq2");
        loadRequiredReal(Params::vq3, vq3_, "vq3");
        loadRequiredReal(Params::lq3, lq3_, "lq3");
        loadRequiredReal(Params::vq4, vq4_, "vq4");
        loadRequiredReal(Params::lq4, lq4_, "lq4");
        loadRequiredReal(Params::vp1, vp1_, "vp1");
        loadRequiredReal(Params::lp1, lp1_, "lp1");
        loadRequiredReal(Params::vp2, vp2_, "vp2");
        loadRequiredReal(Params::lp2, lp2_, "lp2");
        loadRequiredReal(Params::vp3, vp3_, "vp3");
        loadRequiredReal(Params::lp3, lp3_, "lp3");
        loadRequiredReal(Params::vp4, vp4_, "vp4");
        loadRequiredReal(Params::lp4, lp4_, "lp4");
        [[maybe_unused]] bool has_Thld2 = false;
        loadOptionalReal(Params::Thld2, Thld2_, "Thld2", has_Thld2);

        loadOptionalInitReal(External::PE, pe_init_, "pe", has_pe_init_);
        loadOptionalInitReal(External::QGEN, qgen_init_, "qgen", has_qgen_init_);
        loadOptionalInitReal(External::OMEGA, omega_init_, "omega", has_omega_init_);
        loadOptionalInitReal(External::QEXT, qext_init_, "qext", has_qext_init_);
        loadOptionalInitReal(External::PFAREF, pfaref_init_, "pfaref", has_pfaref_init_);
        loadOptionalInitReal(External::PREF, pref_init_, "pref", has_pref_init_);

        if (data.ports.contains(Ports::bus))
        {
          bus_id_ = data.ports.at(Ports::bus);
        }
        else
        {
          Log::error() << "Reeca: missing required port 'bus'\n";
          ++parameter_error_count_;
        }

        setDerivedParameters();
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* Reeca<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void Reeca<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;
        auto index     = [](ReecaInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::iqcmd, [this, index]
                      { return y_[index(ReecaInternalVariables::IQCMD)]; });
        monitor_->set(Variable::ipcmd, [this, index]
                      { return y_[index(ReecaInternalVariables::IPCMD)]; });
        monitor_->set(Variable::vmeas, [this, index]
                      { return y_[index(ReecaInternalVariables::VMEAS)]; });
        monitor_->set(Variable::pmeas, [this, index]
                      { return y_[index(ReecaInternalVariables::PMEAS)]; });
        monitor_->set(Variable::piq, [this, index]
                      { return y_[index(ReecaInternalVariables::XPIQ)]; });
        monitor_->set(Variable::piv, [this, index]
                      { return y_[index(ReecaInternalVariables::XPIV)]; });
        monitor_->set(Variable::qv, [this, index]
                      { return y_[index(ReecaInternalVariables::QV)]; });
        monitor_->set(Variable::pord, [this, index]
                      { return y_[index(ReecaInternalVariables::PORD)]; });
        monitor_->set(Variable::qref, [this, index]
                      { return y_[index(ReecaInternalVariables::QREF)]; });
        monitor_->set(Variable::sdip, [this, index]
                      { return y_[index(ReecaInternalVariables::SDIP)]; });
        monitor_->set(Variable::iqmax, [this, index]
                      { return y_[index(ReecaInternalVariables::IQMAX)]; });
        monitor_->set(Variable::ipmax, [this, index]
                      { return y_[index(ReecaInternalVariables::IPMAX)]; });
        monitor_->set(Variable::iqv, [this, index]
                      { return y_[index(ReecaInternalVariables::IQV)]; });
        monitor_->set(Variable::vqctrl, [this, index]
                      { return y_[index(ReecaInternalVariables::VPIQ)]; });
        monitor_->set(Variable::iqbase, [this, index]
                      { return y_[index(ReecaInternalVariables::IQBASE)]; });
      }

      template <class ScalarT, typename IdxT>
      int Reeca<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Reeca<ScalarT, IdxT>::allocate()
      {
        size_     = static_cast<IdxT>(ReecaInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        f_.assign(size, ScalarT{0});
        y_.assign(size, ScalarT{0});
        yp_.assign(size, ScalarT{0});
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

        static constexpr auto IQCMD = ReecaInternalVariables::IQCMD;
        static constexpr auto IPCMD = ReecaInternalVariables::IPCMD;

        if (signals_.template isAssigned<IQCMD>())
        {
          constexpr auto local_index = static_cast<IdxT>(IQCMD);
          signals_.template getSignalNode<IQCMD>()->set(&y_[static_cast<size_t>(IQCMD)],
                                                        &(this->getVariableIndex(local_index)));
        }

        if (signals_.template isAssigned<IPCMD>())
        {
          constexpr auto local_index = static_cast<IdxT>(IPCMD);
          signals_.template getSignalNode<IPCMD>()->set(&y_[static_cast<size_t>(IPCMD)],
                                                        &(this->getVariableIndex(local_index)));
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Reeca<ScalarT, IdxT>::verify() const
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

        check(Sbase_ > ZERO<RealT>, "Sbase must be positive");
        check(Trv_ >= ZERO<RealT>, "Trv must be non-negative");
        check(Tp_ >= ZERO<RealT>, "Tp must be non-negative");
        check(ZERO<RealT> <= Vdip_ && Vdip_ < Vup_, "Vdip/Vup must satisfy 0 <= Vdip < Vup");
        check(Dbd1_ <= ZERO<RealT> && ZERO<RealT> <= Dbd2_, "Dbd1/Dbd2 must satisfy Dbd1 <= 0 <= Dbd2");
        check(Iql1_ <= Iqh1_, "Iql1 must be less than or equal to Iqh1");
        check(Thld_ == ZERO<RealT>, "Thld must be zero in this REECA implementation");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax");
        check(Tiq_ > ZERO<RealT>, "Tiq must be positive");
        check(Tpord_ > ZERO<RealT>, "Tpord must be positive");
        check(dPmin_ < ZERO<RealT> && ZERO<RealT> < dPmax_, "dPmin < 0 < dPmax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");
        check(Imax_ >= ZERO<RealT>, "Imax must be non-negative");
        check(ZERO<RealT> <= vq1_ && vq1_ < vq2_ && vq2_ < vq3_ && vq3_ < vq4_,
              "vq breakpoints must satisfy 0 <= vq1 < vq2 < vq3 < vq4");
        check(lq1_ >= ZERO<RealT> && lq2_ >= ZERO<RealT> && lq3_ >= ZERO<RealT> && lq4_ >= ZERO<RealT>,
              "lq VDL limits must be non-negative");
        check(ZERO<RealT> <= vp1_ && vp1_ < vp2_ && vp2_ < vp3_ && vp3_ < vp4_,
              "vp breakpoints must satisfy 0 <= vp1 < vp2 < vp3 < vp4");
        check(lp1_ >= ZERO<RealT> && lp2_ >= ZERO<RealT> && lp3_ >= ZERO<RealT> && lp4_ >= ZERO<RealT>,
              "lp VDL limits must be non-negative");
        check(Thld2_ == ZERO<RealT>, "Thld2 must be zero in this REECA implementation");

        static constexpr auto PE     = ReecaExternalVariables::PE;
        static constexpr auto QGEN   = ReecaExternalVariables::QGEN;
        static constexpr auto OMEGA  = ReecaExternalVariables::OMEGA;
        static constexpr auto QEXT   = ReecaExternalVariables::QEXT;
        static constexpr auto PFAREF = ReecaExternalVariables::PFAREF;
        static constexpr auto PREF   = ReecaExternalVariables::PREF;

        if (!isSignalLinked<PE>() && !has_pe_init_)
        {
          Log::error() << "Reeca: pe requires a linked signal or init.pe\n";
          ret += 1;
        }

        if (!isSignalLinked<QGEN>() && !has_qgen_init_)
        {
          Log::error() << "Reeca: qgen requires a linked signal or init.qgen\n";
          ret += 1;
        }

        if (signals_.template isAttached<OMEGA>()
            && !signals_.template isLinked<OMEGA>() && !has_omega_init_)
        {
          Log::error() << "Reeca: omega signal attached with no linked source or init.omega\n";
          ret += 1;
        }

        if (signals_.template isAttached<QEXT>()
            && !signals_.template isLinked<QEXT>() && !has_qext_init_)
        {
          Log::error() << "Reeca: qext signal attached with no linked source or init.qext\n";
          ret += 1;
        }

        if (spf_ && !isSignalLinked<PFAREF>() && !has_pfaref_init_)
        {
          Log::error() << "Reeca: pfaref requires a linked signal or init.pfaref when spf is true\n";
          ret += 1;
        }
        else if (signals_.template isAttached<PFAREF>()
                 && !signals_.template isLinked<PFAREF>() && !has_pfaref_init_)
        {
          Log::error() << "Reeca: pfaref signal attached with no linked source or init.pfaref\n";
          ret += 1;
        }

        if (signals_.template isAttached<PREF>()
            && !signals_.template isLinked<PREF>() && !has_pref_init_)
        {
          Log::error() << "Reeca: pref signal attached with no linked source or init.pref\n";
          ret += 1;
        }

        return ret;
      }

      template <class ScalarT, typename IdxT>
      int Reeca<ScalarT, IdxT>::initialize()
      {
        if (verify() > 0)
        {
          Log::error() << "Reeca: cannot initialize with invalid configuration\n";
          return 1;
        }

        const auto VMEAS      = static_cast<size_t>(ReecaInternalVariables::VMEAS);
        const auto PMEAS      = static_cast<size_t>(ReecaInternalVariables::PMEAS);
        const auto XPIQ       = static_cast<size_t>(ReecaInternalVariables::XPIQ);
        const auto XPIV       = static_cast<size_t>(ReecaInternalVariables::XPIV);
        const auto QV         = static_cast<size_t>(ReecaInternalVariables::QV);
        const auto PORD       = static_cast<size_t>(ReecaInternalVariables::PORD);
        const auto VT         = static_cast<size_t>(ReecaInternalVariables::VT);
        const auto VMEAS_SAFE = static_cast<size_t>(ReecaInternalVariables::VMEAS_SAFE);
        const auto SDIP       = static_cast<size_t>(ReecaInternalVariables::SDIP);
        const auto VERR       = static_cast<size_t>(ReecaInternalVariables::VERR);
        const auto IQV        = static_cast<size_t>(ReecaInternalVariables::IQV);
        const auto QREF       = static_cast<size_t>(ReecaInternalVariables::QREF);
        const auto EQ         = static_cast<size_t>(ReecaInternalVariables::EQ);
        const auto VPIQ       = static_cast<size_t>(ReecaInternalVariables::VPIQ);
        const auto EPIV       = static_cast<size_t>(ReecaInternalVariables::EPIV);
        const auto FPORD      = static_cast<size_t>(ReecaInternalVariables::FPORD);
        const auto RPORD      = static_cast<size_t>(ReecaInternalVariables::RPORD);
        const auto IQCIRC     = static_cast<size_t>(ReecaInternalVariables::IQCIRC);
        const auto IPCIRC     = static_cast<size_t>(ReecaInternalVariables::IPCIRC);
        const auto IQMAX      = static_cast<size_t>(ReecaInternalVariables::IQMAX);
        const auto IPMAX      = static_cast<size_t>(ReecaInternalVariables::IPMAX);
        const auto IQBASE     = static_cast<size_t>(ReecaInternalVariables::IQBASE);
        const auto IQRAW      = static_cast<size_t>(ReecaInternalVariables::IQRAW);
        const auto IQCMD      = static_cast<size_t>(ReecaInternalVariables::IQCMD);
        const auto IPCMD      = static_cast<size_t>(ReecaInternalVariables::IPCMD);

        static constexpr auto PE     = ReecaExternalVariables::PE;
        static constexpr auto QGEN   = ReecaExternalVariables::QGEN;
        static constexpr auto OMEGA  = ReecaExternalVariables::OMEGA;
        static constexpr auto QEXT   = ReecaExternalVariables::QEXT;
        static constexpr auto PFAREF = ReecaExternalVariables::PFAREF;
        static constexpr auto PREF   = ReecaExternalVariables::PREF;

        const ScalarT vr  = Vr();
        const ScalarT vi  = Vi();
        const ScalarT vt0 = std::sqrt(vr * vr + vi * vi);

        if (vt0 <= ZERO<RealT>)
        {
          Log::error() << "Reeca: terminal voltage magnitude must be positive at initialization\n";
          return 1;
        }

        if (!(static_cast<ScalarT>(Vdip_) < vt0 && vt0 < static_cast<ScalarT>(Vup_)))
        {
          Log::error() << "Reeca: standard initialization requires Vdip < Vt0 < Vup\n";
          return 1;
        }

        const ScalarT pe0_system =
            isSignalLinked<PE>() ? signals_.template readExternalVariable<PE>()
                                 : static_cast<ScalarT>(pe_init_);
        const ScalarT qgen0_system =
            isSignalLinked<QGEN>() ? signals_.template readExternalVariable<QGEN>()
                                   : static_cast<ScalarT>(qgen_init_);

        const ScalarT omega0 =
            isSignalLinked<OMEGA>() ? signals_.template readExternalVariable<OMEGA>()
                                    : static_cast<ScalarT>(has_omega_init_ ? omega_init_ : ZERO<RealT>);
        const ScalarT qext0_system =
            isSignalLinked<QEXT>() ? signals_.template readExternalVariable<QEXT>()
                                   : (has_qext_init_ ? static_cast<ScalarT>(qext_init_) : qgen0_system);
        const ScalarT pfaref0 =
            isSignalLinked<PFAREF>() ? signals_.template readExternalVariable<PFAREF>()
                                     : static_cast<ScalarT>(has_pfaref_init_ ? pfaref_init_ : ZERO<RealT>);
        const ScalarT pref_denominator = ONE<RealT> + sP_on_ * omega0;
        if (pref_denominator == ZERO<RealT>)
        {
          Log::error() << "Reeca: pref initialization denominator is zero\n";
          return 1;
        }
        const ScalarT pref0_system =
            isSignalLinked<PREF>() ? signals_.template readExternalVariable<PREF>()
                                   : (has_pref_init_ ? static_cast<ScalarT>(pref_init_)
                                                     : pe0_system / pref_denominator);

        pe_set_     = pe0_system;
        qgen_set_   = qgen0_system;
        omega_set_  = omega0;
        qext_set_   = qext0_system;
        pfaref_set_ = pfaref0;
        pref_set_   = pref0_system;

        if (!has_Vref0_)
        {
          Vref0_ = static_cast<RealT>(vt0);
        }

        const ScalarT pe0    = toReecaBase(pe0_system);
        const ScalarT qgen0  = toReecaBase(qgen0_system);
        const ScalarT qext0  = toReecaBase(qext0_system);
        const ScalarT pref0  = toReecaBase(pref0_system);
        const ScalarT vsafe0 = Math::max(vt0, static_cast<RealT>(0.01));
        const ScalarT sdip0  = Math::outside(vt0, Vdip_, Vup_);
        const ScalarT verr0  = Math::deadband(static_cast<ScalarT>(Vref0_) - vt0, Dbd1_, Dbd2_);
        const ScalarT iqv0   = Math::clamp(Kqv_ * verr0, Iql1_, Iqh1_);
        const ScalarT qref0  = spf_on_ * pe0 * std::tan(pfaref0) + spf_off_ * qext0;

        if (qref0 < static_cast<ScalarT>(Qmin_) || qref0 > static_cast<ScalarT>(Qmax_))
        {
          Log::error() << "Reeca: standard initialization requires Qref0 within Qmin/Qmax\n";
          return 1;
        }

        const ScalarT eq0 = Math::clamp(qref0, Qmin_, Qmax_) - qgen0;
        const RealT   tol = static_cast<RealT>(1.0e-9);
        if (eq0 * eq0 > tol * tol)
        {
          Log::error() << "Reeca: standard initialization requires limited Qref0 equal Qgen0\n";
          return 1;
        }

        const ScalarT vpiq0 =
            sV_ ? vt0 : static_cast<ScalarT>(HALF<RealT> * (Vmin_ + Vmax_));
        if (vpiq0 < static_cast<ScalarT>(Vmin_) || vpiq0 > static_cast<ScalarT>(Vmax_))
        {
          Log::error() << "Reeca: standard initialization requires VPIQ0 within Vmin/Vmax\n";
          return 1;
        }

        const ScalarT epiv0 = sV_on_ * vpiq0 + sV_off_ * (qref0 + Vref1_) - vt0;
        if (epiv0 * epiv0 > tol * tol)
        {
          Log::error() << "Reeca: standard initialization requires ePIV0 equal zero\n";
          return 1;
        }

        const ScalarT pord0 = (ONE<RealT> + sP_on_ * omega0) * pref0;
        if (pord0 < static_cast<ScalarT>(Pmin_) || pord0 > static_cast<ScalarT>(Pmax_))
        {
          Log::error() << "Reeca: standard initialization requires Pord0 within Pmin/Pmax\n";
          return 1;
        }

        const ScalarT fpord0 = ((ONE<RealT> + sP_on_ * omega0) * pref0 - pord0) / Tpord_;
        const ScalarT rpord0 = Math::clamp(fpord0, dPmin_, dPmax_);
        const ScalarT qv0    = qref0 / vsafe0;
        const ScalarT gq0    = vdlQ(vt0);
        const ScalarT gp0    = vdlP(vt0);

        const ScalarT iq_target = qgen0 / vsafe0;
        const ScalarT ip_target = pord0 / vsafe0;

        ScalarT iqcirc0{0};
        ScalarT ipcirc0{0};
        ScalarT iqmax0{0};
        ScalarT ipmax0{0};
        ScalarT iqbase0{iq_target};
        ScalarT iqraw0{0};
        ScalarT iqcmd0{0};
        ScalarT ipcmd0{0};

        if (sPQ_)
        {
          ipcirc0 = Imax_;
          ipmax0  = Math::min(gp0, ipcirc0);
          if (ip_target < ZERO<RealT> || ip_target > ipmax0)
          {
            Log::error() << "Reeca: standard initialization requires Ip0 within active-current limits\n";
            return 1;
          }
          ipcmd0 = Math::clamp(ip_target, ZERO<RealT>, ipmax0);

          const ScalarT iq_radicand = static_cast<ScalarT>(Imax_ * Imax_) - ipcmd0 * ipcmd0;
          if (iq_radicand < ZERO<RealT>)
          {
            Log::error() << "Reeca: standard initialization has negative reactive-current circle radicand\n";
            return 1;
          }
          iqcirc0 = std::sqrt(iq_radicand);
          iqmax0  = Math::min(gq0, iqcirc0);
          if (iq_target < -iqmax0 || iq_target > iqmax0)
          {
            Log::error() << "Reeca: standard initialization requires Iq0 within reactive-current limits\n";
            return 1;
          }
          iqraw0 = sQ_on_ * iqbase0 + sQ_off_ * qv0 + sdip0 * iqv0;
          if (iqraw0 < -iqmax0 || iqraw0 > iqmax0)
          {
            Log::error() << "Reeca: standard initialization requires raw Iq0 within reactive-current limits\n";
            return 1;
          }
          iqcmd0 = Math::clamp(iqraw0, -iqmax0, iqmax0);
        }
        else
        {
          iqcirc0 = Imax_;
          iqmax0  = Math::min(gq0, iqcirc0);
          if (iq_target < -iqmax0 || iq_target > iqmax0)
          {
            Log::error() << "Reeca: standard initialization requires Iq0 within reactive-current limits\n";
            return 1;
          }
          iqraw0 = sQ_on_ * iqbase0 + sQ_off_ * qv0 + sdip0 * iqv0;
          if (iqraw0 < -iqmax0 || iqraw0 > iqmax0)
          {
            Log::error() << "Reeca: standard initialization requires raw Iq0 within reactive-current limits\n";
            return 1;
          }
          iqcmd0 = Math::clamp(iqraw0, -iqmax0, iqmax0);

          const ScalarT ip_radicand = static_cast<ScalarT>(Imax_ * Imax_) - iqcmd0 * iqcmd0;
          if (ip_radicand < ZERO<RealT>)
          {
            Log::error() << "Reeca: standard initialization has negative active-current circle radicand\n";
            return 1;
          }
          ipcirc0 = std::sqrt(ip_radicand);
          ipmax0  = Math::min(gp0, ipcirc0);
          if (ip_target < ZERO<RealT> || ip_target > ipmax0)
          {
            Log::error() << "Reeca: standard initialization requires Ip0 within active-current limits\n";
            return 1;
          }
          ipcmd0 = Math::clamp(ip_target, ZERO<RealT>, ipmax0);
        }

        std::fill(y_.begin(), y_.end(), ZERO<RealT>);
        std::fill(yp_.begin(), yp_.end(), ZERO<RealT>);

        y_[VMEAS]      = vt0;
        y_[PMEAS]      = pe0;
        y_[XPIQ]       = vpiq0 - Kqp_ * eq0;
        y_[XPIV]       = iqbase0 - Kvp_ * epiv0;
        y_[QV]         = qv0;
        y_[PORD]       = pord0;
        y_[VT]         = vt0;
        y_[VMEAS_SAFE] = vsafe0;
        y_[SDIP]       = sdip0;
        y_[VERR]       = verr0;
        y_[IQV]        = iqv0;
        y_[QREF]       = qref0;
        y_[EQ]         = eq0;
        y_[VPIQ]       = vpiq0;
        y_[EPIV]       = epiv0;
        y_[FPORD]      = fpord0;
        y_[RPORD]      = rpord0;
        y_[IQCIRC]     = iqcirc0;
        y_[IPCIRC]     = ipcirc0;
        y_[IQMAX]      = iqmax0;
        y_[IPMAX]      = ipmax0;
        y_[IQBASE]     = iqbase0;
        y_[IQRAW]      = iqraw0;
        y_[IQCMD]      = iqcmd0;
        y_[IPCMD]      = ipcmd0;

        if (!ws_.empty())
        {
          ws_[static_cast<size_t>(PE)]     = pe0_system;
          ws_[static_cast<size_t>(QGEN)]   = qgen0_system;
          ws_[static_cast<size_t>(OMEGA)]  = omega_set_;
          ws_[static_cast<size_t>(QEXT)]   = qext_set_;
          ws_[static_cast<size_t>(PFAREF)] = pfaref_set_;
          ws_[static_cast<size_t>(PREF)]   = pref_set_;
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Reeca<ScalarT, IdxT>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(ReecaInternalVariables::VMEAS)] = (Trv_ > ZERO<RealT>);
        tag_[static_cast<size_t>(ReecaInternalVariables::PMEAS)] = (Tp_ > ZERO<RealT>);
        tag_[static_cast<size_t>(ReecaInternalVariables::XPIQ)]  = true;
        tag_[static_cast<size_t>(ReecaInternalVariables::XPIV)]  = true;
        tag_[static_cast<size_t>(ReecaInternalVariables::QV)]    = true;
        tag_[static_cast<size_t>(ReecaInternalVariables::PORD)]  = true;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Reeca<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT* y,
          ScalarT* yp,
          ScalarT* wb,
          ScalarT* ws,
          ScalarT* f)
      {
        const auto VMEAS      = static_cast<size_t>(ReecaInternalVariables::VMEAS);
        const auto PMEAS      = static_cast<size_t>(ReecaInternalVariables::PMEAS);
        const auto XPIQ       = static_cast<size_t>(ReecaInternalVariables::XPIQ);
        const auto XPIV       = static_cast<size_t>(ReecaInternalVariables::XPIV);
        const auto QV         = static_cast<size_t>(ReecaInternalVariables::QV);
        const auto PORD       = static_cast<size_t>(ReecaInternalVariables::PORD);
        const auto VT         = static_cast<size_t>(ReecaInternalVariables::VT);
        const auto VMEAS_SAFE = static_cast<size_t>(ReecaInternalVariables::VMEAS_SAFE);
        const auto SDIP       = static_cast<size_t>(ReecaInternalVariables::SDIP);
        const auto VERR       = static_cast<size_t>(ReecaInternalVariables::VERR);
        const auto IQV        = static_cast<size_t>(ReecaInternalVariables::IQV);
        const auto QREF       = static_cast<size_t>(ReecaInternalVariables::QREF);
        const auto EQ         = static_cast<size_t>(ReecaInternalVariables::EQ);
        const auto VPIQ       = static_cast<size_t>(ReecaInternalVariables::VPIQ);
        const auto EPIV       = static_cast<size_t>(ReecaInternalVariables::EPIV);
        const auto FPORD      = static_cast<size_t>(ReecaInternalVariables::FPORD);
        const auto RPORD      = static_cast<size_t>(ReecaInternalVariables::RPORD);
        const auto IQCIRC     = static_cast<size_t>(ReecaInternalVariables::IQCIRC);
        const auto IPCIRC     = static_cast<size_t>(ReecaInternalVariables::IPCIRC);
        const auto IQMAX      = static_cast<size_t>(ReecaInternalVariables::IQMAX);
        const auto IPMAX      = static_cast<size_t>(ReecaInternalVariables::IPMAX);
        const auto IQBASE     = static_cast<size_t>(ReecaInternalVariables::IQBASE);
        const auto IQRAW      = static_cast<size_t>(ReecaInternalVariables::IQRAW);
        const auto IQCMD      = static_cast<size_t>(ReecaInternalVariables::IQCMD);
        const auto IPCMD      = static_cast<size_t>(ReecaInternalVariables::IPCMD);

        const auto PE     = static_cast<size_t>(ReecaExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecaExternalVariables::QGEN);
        const auto OMEGA  = static_cast<size_t>(ReecaExternalVariables::OMEGA);
        const auto QEXT   = static_cast<size_t>(ReecaExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecaExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecaExternalVariables::PREF);

        const ScalarT vr     = wb[0];
        const ScalarT vi     = wb[1];
        const ScalarT vt     = y[VT];
        const ScalarT pe     = toReecaBase(ws[PE]);
        const ScalarT qgen   = toReecaBase(ws[QGEN]);
        const ScalarT qext   = toReecaBase(ws[QEXT]);
        const ScalarT pref   = toReecaBase(ws[PREF]);
        const ScalarT omega  = ws[OMEGA];
        const ScalarT pfaref = ws[PFAREF];

        const ScalarT vsafe = Math::max(y[VMEAS], static_cast<RealT>(0.01));
        const ScalarT sdip  = Math::outside(vt, Vdip_, Vup_);
        const ScalarT verr  = Math::deadband(static_cast<ScalarT>(Vref0_) - y[VMEAS], Dbd1_, Dbd2_);
        const ScalarT iqv   = Math::clamp(Kqv_ * verr, Iql1_, Iqh1_);
        const ScalarT qref  = spf_on_ * y[PMEAS] * std::tan(pfaref) + spf_off_ * qext;
        const ScalarT eq    = Math::clamp(qref, Qmin_, Qmax_) - qgen;
        const ScalarT vpiq  = Math::clamp(Kqp_ * eq + y[XPIQ], Vmin_, Vmax_);
        const ScalarT epiv  = sV_on_ * vpiq + sV_off_ * (qref + Vref1_) - y[VMEAS];
        const ScalarT fpord = ((ONE<RealT> + sP_on_ * omega) * pref - y[PORD]) / Tpord_;
        const ScalarT rpord = Math::clamp(fpord, dPmin_, dPmax_);

        const ScalarT gq = static_cast<ScalarT>(lq1_)
                           + Math::linseg(y[VMEAS], vq1_, vq2_, lq2_ - lq1_)
                           + Math::linseg(y[VMEAS], vq2_, vq3_, lq3_ - lq2_)
                           + Math::linseg(y[VMEAS], vq3_, vq4_, lq4_ - lq3_);
        const ScalarT gp = static_cast<ScalarT>(lp1_)
                           + Math::linseg(y[VMEAS], vp1_, vp2_, lp2_ - lp1_)
                           + Math::linseg(y[VMEAS], vp2_, vp3_, lp3_ - lp2_)
                           + Math::linseg(y[VMEAS], vp3_, vp4_, lp4_ - lp3_);

        const ScalarT xpiv_rate      = Kvi_ * epiv;
        const ScalarT xpiv_above_min = Math::sigmoid(y[IQBASE] + y[IQMAX]);
        const ScalarT xpiv_below_max = Math::sigmoid(y[IQMAX] - y[IQBASE]);
        const ScalarT xpiv_aw =
            (xpiv_above_min * xpiv_below_max
             + (ONE<RealT> - xpiv_below_max) * Math::sigmoid(-xpiv_rate)
             + (ONE<RealT> - xpiv_above_min) * Math::sigmoid(xpiv_rate))
            * xpiv_rate;

        const ScalarT iqbase_pre = Kvp_ * y[EPIV] + y[XPIV];
        const ScalarT iqbase =
            -y[IQMAX] + Math::ramp(iqbase_pre + y[IQMAX]) - Math::ramp(iqbase_pre - y[IQMAX]);
        const ScalarT iqcmd =
            -y[IQMAX] + Math::ramp(y[IQRAW] + y[IQMAX]) - Math::ramp(y[IQRAW] - y[IQMAX]);
        const ScalarT ipcmd_pre = y[PORD] / y[VMEAS_SAFE];
        const ScalarT ipcmd     = Math::ramp(ipcmd_pre) - Math::ramp(ipcmd_pre - y[IPMAX]);

        f[VMEAS]      = -Trv_ * yp[VMEAS] - y[VMEAS] + vt;
        f[PMEAS]      = -Tp_ * yp[PMEAS] - y[PMEAS] + pe;
        f[XPIQ]       = -yp[XPIQ] + (ONE<RealT> - sdip) * Math::antiwindup(vpiq, Kqi_ * eq, Vmin_, Vmax_);
        f[XPIV]       = -yp[XPIV] + (ONE<RealT> - sdip) * xpiv_aw;
        f[QV]         = -Tiq_ * yp[QV] - (ONE<RealT> - sdip) * y[QV] + (ONE<RealT> - sdip) * qref / vsafe;
        f[PORD]       = -yp[PORD] + (ONE<RealT> - sdip) * Math::antiwindup(y[PORD], rpord, Pmin_, Pmax_);
        f[VT]         = -vt * vt + vr * vr + vi * vi;
        f[VMEAS_SAFE] = -y[VMEAS_SAFE] + vsafe;
        f[SDIP]       = -y[SDIP] + sdip;
        f[VERR]       = -y[VERR] + verr;
        f[IQV]        = -y[IQV] + iqv;
        f[QREF]       = -y[QREF] + qref;
        f[EQ]         = -y[EQ] + eq;
        f[VPIQ]       = -y[VPIQ] + vpiq;
        f[EPIV]       = -y[EPIV] + epiv;
        f[FPORD]      = -y[FPORD] + fpord;
        f[RPORD]      = -y[RPORD] + rpord;
        f[IQCIRC]     = -y[IQCIRC] * y[IQCIRC] + Imax_ * Imax_ - sPQ_on_ * y[IPCMD] * y[IPCMD];
        f[IPCIRC]     = -y[IPCIRC] * y[IPCIRC] + Imax_ * Imax_ - sPQ_off_ * y[IQCMD] * y[IQCMD];
        f[IQMAX]      = -y[IQMAX] + Math::min(gq, y[IQCIRC]);
        f[IPMAX]      = -y[IPMAX] + Math::min(gp, y[IPCIRC]);
        f[IQBASE]     = -y[IQBASE] + iqbase;
        f[IQRAW]      = -y[IQRAW] + sQ_on_ * y[IQBASE] + sQ_off_ * y[QV] + y[SDIP] * y[IQV];
        f[IQCMD]      = -y[IQCMD] + iqcmd;
        f[IPCMD]      = -y[IPCMD] + ipcmd;

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Reeca<ScalarT, IdxT>::evaluateResidual()
      {
        static constexpr auto PE     = ReecaExternalVariables::PE;
        static constexpr auto QGEN   = ReecaExternalVariables::QGEN;
        static constexpr auto OMEGA  = ReecaExternalVariables::OMEGA;
        static constexpr auto QEXT   = ReecaExternalVariables::QEXT;
        static constexpr auto PFAREF = ReecaExternalVariables::PFAREF;
        static constexpr auto PREF   = ReecaExternalVariables::PREF;

        std::fill(ws_.begin(), ws_.end(), ZERO<RealT>);
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        ws_[static_cast<size_t>(PE)]     = pe_set_;
        ws_[static_cast<size_t>(QGEN)]   = qgen_set_;
        ws_[static_cast<size_t>(OMEGA)]  = omega_set_;
        ws_[static_cast<size_t>(QEXT)]   = qext_set_;
        ws_[static_cast<size_t>(PFAREF)] = pfaref_set_;
        ws_[static_cast<size_t>(PREF)]   = pref_set_;

        if (signals_.template isAttached<PE>() && signals_.template isLinked<PE>())
        {
          ws_[static_cast<size_t>(PE)] = signals_.template readExternalVariable<PE>();
          ws_indices_[static_cast<size_t>(PE)] =
              signals_.template readExternalVariableIndex<PE>();
        }

        if (signals_.template isAttached<QGEN>() && signals_.template isLinked<QGEN>())
        {
          ws_[static_cast<size_t>(QGEN)] = signals_.template readExternalVariable<QGEN>();
          ws_indices_[static_cast<size_t>(QGEN)] =
              signals_.template readExternalVariableIndex<QGEN>();
        }

        if (signals_.template isAttached<OMEGA>() && signals_.template isLinked<OMEGA>())
        {
          ws_[static_cast<size_t>(OMEGA)] = signals_.template readExternalVariable<OMEGA>();
          ws_indices_[static_cast<size_t>(OMEGA)] =
              signals_.template readExternalVariableIndex<OMEGA>();
        }

        if (signals_.template isAttached<QEXT>() && signals_.template isLinked<QEXT>())
        {
          ws_[static_cast<size_t>(QEXT)] = signals_.template readExternalVariable<QEXT>();
          ws_indices_[static_cast<size_t>(QEXT)] =
              signals_.template readExternalVariableIndex<QEXT>();
        }

        if (signals_.template isAttached<PFAREF>() && signals_.template isLinked<PFAREF>())
        {
          ws_[static_cast<size_t>(PFAREF)] = signals_.template readExternalVariable<PFAREF>();
          ws_indices_[static_cast<size_t>(PFAREF)] =
              signals_.template readExternalVariableIndex<PFAREF>();
        }

        if (signals_.template isAttached<PREF>() && signals_.template isLinked<PREF>())
        {
          ws_[static_cast<size_t>(PREF)] = signals_.template readExternalVariable<PREF>();
          ws_indices_[static_cast<size_t>(PREF)] =
              signals_.template readExternalVariableIndex<PREF>();
        }

        wb_[0] = Vr();
        wb_[1] = Vi();

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

        return 0;
      }
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

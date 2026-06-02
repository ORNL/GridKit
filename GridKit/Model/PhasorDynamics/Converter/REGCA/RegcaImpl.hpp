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

      template <class ScalarT, typename IdxT>
      Regca<ScalarT, IdxT>::Regca(bus_type* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Regca<ScalarT, IdxT>::Regca(bus_type* bus, const model_data_type& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Regca<ScalarT, IdxT>::~Regca()
      {
      }

      template <class ScalarT, typename IdxT>
      void Regca<ScalarT, IdxT>::setDerivedParameters()
      {
        Mp_          = static_cast<RealT>(100.0) * Rpmax_;
        use_lvpl_    = ZERO<RealT>;
        bypass_lvpl_ = ONE<RealT>;
        if (sL_)
        {
          use_lvpl_    = ONE<RealT>;
          bypass_lvpl_ = ZERO<RealT>;
        }
        iq_use_upper_      = ZERO<RealT>;
        iq_use_lower_      = ONE<RealT>;
        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);
      }

      template <class ScalarT, typename IdxT>
      ScalarT Regca<ScalarT, IdxT>::activeCurrentLowerRateBound(ScalarT ip) const
      {
        return -Rpmax_ - (Mp_ - Rpmax_) * Math::sigmoid(ip);
      }

      template <class ScalarT, typename IdxT>
      ScalarT Regca<ScalarT, IdxT>::activeCurrentUpperRateBound(ScalarT ip, ScalarT il) const
      {
        const ScalarT sigma_ip = Math::sigmoid(ip);
        return Mp_ * (ONE<RealT> - sigma_ip)
               + Rpmax_ * sigma_ip * (bypass_lvpl_ + use_lvpl_ * Math::sigmoid(il - ip));
      }

      template <class ScalarT, typename IdxT>
      void Regca<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
      {
        using Params = typename model_data_type::Parameters;
        using Ports  = typename model_data_type::Ports;

        parameter_error_count_ = 0;

        auto loadRequiredReal = [&](auto key, RealT& target, const char* name)
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

        auto loadRequiredSwitch = [&](auto key, bool& target, const char* name)
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

        loadRequiredReal(Params::P0, P0_, "P0");
        loadRequiredReal(Params::Q0, Q0_, "Q0");
        loadRequiredReal(Params::mva, mva_base_, "mva");
        loadRequiredReal(Params::Tg, Tg_, "Tg");
        loadRequiredReal(Params::TM, TM_, "TM");
        loadRequiredReal(Params::Rqmax, Rqmax_, "Rqmax");
        loadRequiredReal(Params::Rqmin, Rqmin_, "Rqmin");
        loadRequiredReal(Params::Rpmax, Rpmax_, "Rpmax");
        loadRequiredSwitch(Params::sL, sL_, "sL");
        loadRequiredReal(Params::IL1, IL1_, "IL1");
        loadRequiredReal(Params::VL0, VL0_, "VL0");
        loadRequiredReal(Params::VL1, VL1_, "VL1");
        loadRequiredReal(Params::VA0, VA0_, "VA0");
        loadRequiredReal(Params::VA1, VA1_, "VA1");
        loadRequiredReal(Params::Vhvmax, Vhvmax_, "Vhvmax");

        if (data.ports.contains(Ports::bus))
        {
          bus_id_ = data.ports.at(Ports::bus);
        }

        setDerivedParameters();
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* Regca<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void Regca<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;
        auto index     = [](RegcaInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::ir, [this, index]
                      { return y_[index(RegcaInternalVariables::IR)]; });
        monitor_->set(Variable::ii, [this, index]
                      { return y_[index(RegcaInternalVariables::II)]; });
        monitor_->set(Variable::p, [this, index]
                      { return y_[index(RegcaInternalVariables::PBR)]; });
        monitor_->set(Variable::q, [this, index]
                      { return y_[index(RegcaInternalVariables::QBR)]; });
        monitor_->set(Variable::vt, [this, index]
                      { return y_[index(RegcaInternalVariables::VT)]; });
        monitor_->set(Variable::vm, [this, index]
                      { return y_[index(RegcaInternalVariables::VM)]; });
        monitor_->set(Variable::ip, [this, index]
                      { return y_[index(RegcaInternalVariables::IP)]; });
        monitor_->set(Variable::iq, [this, index]
                      { return y_[index(RegcaInternalVariables::IQ)]; });
        monitor_->set(Variable::iqextra, [this, index]
                      { return y_[index(RegcaInternalVariables::IQEXTRA)]; });
        monitor_->set(Variable::il, [this, index]
                      { return y_[index(RegcaInternalVariables::IL)]; });
        monitor_->set(Variable::lp, [this, index]
                      { return y_[index(RegcaInternalVariables::LP)]; });
        monitor_->set(Variable::up, [this, index]
                      { return y_[index(RegcaInternalVariables::UP)]; });
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::allocate()
      {
        size_     = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        f_.assign(size, ScalarT{0});
        y_.assign(size, ScalarT{0});
        yp_.assign(size, ScalarT{0});
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

        if (signals_.template isAssigned<RegcaInternalVariables::IR>())
        {
          signals_.template getSignalNode<RegcaInternalVariables::IR>()->set(
              &y_[static_cast<size_t>(RegcaInternalVariables::IR)],
              &(this->getVariableIndex(static_cast<IdxT>(RegcaInternalVariables::IR))));
        }

        if (signals_.template isAssigned<RegcaInternalVariables::II>())
        {
          signals_.template getSignalNode<RegcaInternalVariables::II>()->set(
              &y_[static_cast<size_t>(RegcaInternalVariables::II)],
              &(this->getVariableIndex(static_cast<IdxT>(RegcaInternalVariables::II))));
        }

        if (signals_.template isAssigned<RegcaInternalVariables::PBR>())
        {
          signals_.template getSignalNode<RegcaInternalVariables::PBR>()->set(
              &y_[static_cast<size_t>(RegcaInternalVariables::PBR)],
              &(this->getVariableIndex(static_cast<IdxT>(RegcaInternalVariables::PBR))));
        }

        if (signals_.template isAssigned<RegcaInternalVariables::QBR>())
        {
          signals_.template getSignalNode<RegcaInternalVariables::QBR>()->set(
              &y_[static_cast<size_t>(RegcaInternalVariables::QBR)],
              &(this->getVariableIndex(static_cast<IdxT>(RegcaInternalVariables::QBR))));
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::verify() const
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
        check(Tg_ > ZERO<RealT>, "Tg must be positive");
        check(TM_ > ZERO<RealT>, "TM must be positive");
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

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::initialize()
      {
        if (bus_ == nullptr)
        {
          Log::error() << "Regca: cannot initialize with null bus\n";
          return 1;
        }

        if (parameter_error_count_ > 0 || mva_base_ <= ZERO<RealT> || Tg_ <= ZERO<RealT>
            || TM_ <= ZERO<RealT> || Rpmax_ <= ZERO<RealT>
            || !(Rqmin_ < ZERO<RealT> && ZERO<RealT> < Rqmax_)
            || IL1_ < ZERO<RealT> || !(ZERO<RealT> <= VL0_ && VL0_ < VL1_)
            || !(ZERO<RealT> <= VA0_ && VA0_ < VA1_) || Vhvmax_ <= ZERO<RealT>)
        {
          Log::error() << "Regca: cannot initialize with invalid parameters\n";
          return 1;
        }

        const auto VM      = static_cast<size_t>(RegcaInternalVariables::VM);
        const auto IQ      = static_cast<size_t>(RegcaInternalVariables::IQ);
        const auto IP      = static_cast<size_t>(RegcaInternalVariables::IP);
        const auto VT      = static_cast<size_t>(RegcaInternalVariables::VT);
        const auto II      = static_cast<size_t>(RegcaInternalVariables::II);
        const auto IQEXTRA = static_cast<size_t>(RegcaInternalVariables::IQEXTRA);
        const auto IL      = static_cast<size_t>(RegcaInternalVariables::IL);
        const auto IR      = static_cast<size_t>(RegcaInternalVariables::IR);
        const auto LP      = static_cast<size_t>(RegcaInternalVariables::LP);
        const auto UP      = static_cast<size_t>(RegcaInternalVariables::UP);
        const auto PBR     = static_cast<size_t>(RegcaInternalVariables::PBR);
        const auto QBR     = static_cast<size_t>(RegcaInternalVariables::QBR);

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();
        const ScalarT vt = std::sqrt(vr * vr + vi * vi);

        if (vt <= ZERO<RealT>)
        {
          Log::error() << "Regca: terminal voltage magnitude must be positive at initialization\n";
          return 1;
        }
        if (vt >= Vhvmax_)
        {
          Log::error() << "Regca: terminal voltage magnitude must be below Vhvmax at initialization\n";
          return 1;
        }

        // REGCA owns the network terminal and establishes the initial
        // converter operating point from the power-flow injection. Controller
        // command ports may not have been initialized yet, so initialization
        // resolves the commands from P0/Q0 and publishes them to attached
        // ports below. P0/Q0 are given on the system base; toComponentBase
        // converts them to the converter base the internal states use.
        const ScalarT lvacm = Math::linseg(vt, VA0_, VA1_, ONE<RealT>);

        if (P0_ != ZERO<RealT> && (vt <= VA0_ || lvacm <= ZERO<RealT>) )
        {
          Log::error() << "Regca: LVACM gain is zero with nonzero initial active power\n";
          return 1;
        }

        ScalarT ipcmd0{ZERO<RealT>};
        if (P0_ != ZERO<RealT>)
        {
          ipcmd0 = toComponentBase(static_cast<ScalarT>(P0_) / vt) / lvacm;
        }
        const ScalarT iqcmd0 = toComponentBase(static_cast<ScalarT>(Q0_) / vt);

        const ScalarT iqextra0{ZERO<RealT>};
        const ScalarT qnet0 = iqcmd0 - iqextra0;

        y_[VM]      = vt;
        y_[VT]      = vt;
        y_[IP]      = ipcmd0;
        y_[IQ]      = iqcmd0;
        y_[IQEXTRA] = iqextra0;
        y_[IL]      = Math::linseg(vt, VL0_, VL1_, IL1_);
        y_[IR]      = (qnet0 * vi + lvacm * ipcmd0 * vr) / vt;
        y_[II]      = (-qnet0 * vr + lvacm * ipcmd0 * vi) / vt;
        y_[LP]      = activeCurrentLowerRateBound(y_[IP]);
        y_[UP]      = activeCurrentUpperRateBound(y_[IP], y_[IL]);
        y_[PBR]     = toSystemBase(vr * y_[IR] + vi * y_[II]);
        y_[QBR]     = toSystemBase(vi * y_[IR] - vr * y_[II]);

        // Retain the resolved commands as the constant source used during
        // residual evaluation when no controller drives the command ports, and
        // select the reactive-current rate-limit branch from the command sign.
        ipcmd_set_    = ipcmd0;
        iqcmd_set_    = iqcmd0;
        iq_use_upper_ = ZERO<RealT>;
        iq_use_lower_ = ONE<RealT>;
        if (static_cast<RealT>(iqcmd0) > ZERO<RealT>)
        {
          iq_use_upper_ = ONE<RealT>;
          iq_use_lower_ = ZERO<RealT>;
        }

        // Seed attached command nodes with the steady-state values. Controller
        // initialization can use these signal values, and unattached ports fall
        // back to the constants stored above.
        if (signals_.template isAttached<RegcaExternalVariables::IPCMD>())
        {
          signals_.template writeExternalVariable<RegcaExternalVariables::IPCMD>(ipcmd0);
        }
        if (signals_.template isAttached<RegcaExternalVariables::IQCMD>())
        {
          signals_.template writeExternalVariable<RegcaExternalVariables::IQCMD>(iqcmd0);
        }

        std::fill(yp_.begin(), yp_.end(), ZERO<RealT>);
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(RegcaInternalVariables::VM)] = true;
        tag_[static_cast<size_t>(RegcaInternalVariables::IQ)] = true;
        tag_[static_cast<size_t>(RegcaInternalVariables::IP)] = true;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Regca<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT* y,
          ScalarT* yp,
          ScalarT* wb,
          ScalarT* ws,
          ScalarT* f)
      {
        const auto VM      = static_cast<size_t>(RegcaInternalVariables::VM);
        const auto IQ      = static_cast<size_t>(RegcaInternalVariables::IQ);
        const auto IP      = static_cast<size_t>(RegcaInternalVariables::IP);
        const auto VT      = static_cast<size_t>(RegcaInternalVariables::VT);
        const auto II      = static_cast<size_t>(RegcaInternalVariables::II);
        const auto IQEXTRA = static_cast<size_t>(RegcaInternalVariables::IQEXTRA);
        const auto IL      = static_cast<size_t>(RegcaInternalVariables::IL);
        const auto IR      = static_cast<size_t>(RegcaInternalVariables::IR);
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
        const ScalarT ii      = y[II];
        const ScalarT iqextra = y[IQEXTRA];
        const ScalarT il      = y[IL];
        const ScalarT ir      = y[IR];
        const ScalarT lp      = y[LP];
        const ScalarT up      = y[UP];
        const ScalarT pbr     = y[PBR];
        const ScalarT qbr     = y[QBR];

        const ScalarT vm_dot = yp[VM];
        const ScalarT iq_dot = yp[IQ];
        const ScalarT ip_dot = yp[IP];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT ipcmd = ws[IPCMD];
        const ScalarT iqcmd = ws[IQCMD];

        const ScalarT iq_error = iqcmd - iq;
        const ScalarT ip_error = ipcmd - ip;

        const ScalarT iq_step =
            iq_use_upper_ * (iq_error - Math::ramp(iq_error - Tg_ * Rqmax_))
            + iq_use_lower_ * (iq_error + Math::ramp(Tg_ * Rqmin_ - iq_error));

        const ScalarT ip_step =
            Tg_ * lp + Math::ramp(ip_error - Tg_ * lp) - Math::ramp(ip_error - Tg_ * up);

        f[VM]               = -TM_ * vm_dot - vm + vt;
        f[IQ]               = -Tg_ * iq_dot + iq_step;
        f[IP]               = -Tg_ * ip_dot + ip_step;
        const ScalarT lvacm = Math::linseg(vt, VA0_, VA1_, ONE<RealT>);
        const ScalarT qnet  = iq - iqextra;

        f[VT]      = -vt * vt + vr * vr + vi * vi;
        f[II]      = -vt * ii - qnet * vr + lvacm * ip * vi;
        f[IQEXTRA] = -iqextra + Math::ramp(vt - Vhvmax_);
        f[IL]      = -il + Math::linseg(vm, VL0_, VL1_, IL1_);
        f[IR]      = -vt * ir + qnet * vi + lvacm * ip * vr;
        f[LP]      = -lp + activeCurrentLowerRateBound(ip);
        f[UP]      = -up + activeCurrentUpperRateBound(ip, il);
        f[PBR]     = -pbr + toSystemBase(vr * ir + vi * ii);
        f[QBR]     = -qbr + toSystemBase(vi * ir - vr * ii);

        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Regca<ScalarT, IdxT>::evaluateBusResidual(
          ScalarT*                  y,
          [[maybe_unused]] ScalarT* yp,
          [[maybe_unused]] ScalarT* wb,
          ScalarT*                  h)
      {
        const auto II = static_cast<size_t>(RegcaInternalVariables::II);
        const auto IR = static_cast<size_t>(RegcaInternalVariables::IR);

        h[0] = toSystemBase(y[IR]);
        h[1] = toSystemBase(y[II]);
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::evaluateResidual()
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

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
        evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());

        Ir() += h_[0];
        Ii() += h_[1];

        return 0;
      }
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit

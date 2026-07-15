/**
 * @file Esdc1aImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the ESDC1A exciter model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1a.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      using Log = ::GridKit::Utilities::Logger;

      template <class ScalarT, typename IdxT>
      Esdc1a<ScalarT, IdxT>::Esdc1a(bus_type* bus)
        : bus_(bus)
      {
        setDerivedParams();
        size_ = static_cast<IdxT>(Esdc1aInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Esdc1a<ScalarT, IdxT>::Esdc1a(bus_type* bus, const model_data_type& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        setDerivedParams();
        initializeMonitor();
        size_ = static_cast<IdxT>(Esdc1aInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Esdc1a<ScalarT, IdxT>::~Esdc1a()
      {
      }

      template <class ScalarT, typename IdxT>
      void Esdc1a<ScalarT, IdxT>::initModelParams(const model_data_type& data)
      {
        using Params = typename model_data_type::Parameters;

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
            Log::error() << "Esdc1a: parameter '" << name << "' must be numeric\n";
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
            target = *bool_value ? ONE<RealT> : ZERO<RealT>;
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
            Log::error() << "Esdc1a: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        auto load_selector = [&](auto key, IdxT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            target = *index_value;
          }
          else if (const auto* real_value = std::get_if<RealT>(&value))
          {
            const RealT rounded = std::round(*real_value);
            if (*real_value >= ZERO<RealT> && *real_value == rounded)
            {
              target = static_cast<IdxT>(rounded);
            }
            else
            {
              Log::error() << "Esdc1a: parameter '" << name << "' must be an integer selector\n";
              ++parameter_error_count_;
            }
          }
          else
          {
            Log::error() << "Esdc1a: parameter '" << name << "' must be an integer selector\n";
            ++parameter_error_count_;
          }
        };

        load_real(Params::Tr, Tr_, "Tr");
        load_real(Params::Ka, Ka_, "Ka");
        load_real(Params::Ta, Ta_, "Ta");
        load_real(Params::Tb, Tb_, "Tb");
        load_real(Params::Tc, Tc_, "Tc");
        load_real(Params::Vrmax, Vrmax_, "Vrmax");
        load_real(Params::Vrmin, Vrmin_, "Vrmin");
        load_real(Params::Ke, Ke_, "Ke");
        load_real(Params::Te, Te_, "Te");
        load_real(Params::Kf, Kf_, "Kf");
        load_real(Params::Tf1, Tf1_, "Tf1");
        load_switch(Params::Spdmlt, spdmlt_, "Spdmlt");
        load_real(Params::E1, E1_, "E1");
        load_real(Params::Se1, Se1_, "Se1");
        load_real(Params::E2, E2_, "E2");
        load_real(Params::Se2, Se2_, "Se2");
        load_selector(Params::UEL, UEL_, "UEL");
        load_switch(Params::exclim, exclim_, "exclim");
      }

      template <class ScalarT, typename IdxT>
      void Esdc1a<ScalarT, IdxT>::setDerivedParams()
      {
        Tr_  = std::max(Tr_, TIME_CONSTANT_MINIMUM);
        Tb_  = std::max(Tb_, TIME_CONSTANT_MINIMUM);
        Tf1_ = std::max(Tf1_, TIME_CONSTANT_MINIMUM);

        sUEL_     = UEL_ >= static_cast<IdxT>(2) ? ONE<RealT> : ZERO<RealT>;
        sUELoff_  = ONE<RealT> - sUEL_;
        slim_     = exclim_;
        slim_off_ = ONE<RealT> - slim_;

        if (Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>)
        {
          SA_ = ZERO<RealT>;
          SB_ = ZERO<RealT>;
          return;
        }
        if (E1_ <= ZERO<RealT> || E2_ <= ZERO<RealT> || E1_ == E2_
            || Se1_ <= ZERO<RealT> || Se2_ <= ZERO<RealT> || Se1_ == Se2_)
        {
          SA_ = ZERO<RealT>;
          SB_ = ZERO<RealT>;
          return;
        }

        const RealT C = std::sqrt(Se2_ / Se1_);
        SA_           = (C * E1_ - E2_) / (C - ONE<RealT>);
        SB_           = Se1_ / ((E1_ - SA_) * (E1_ - SA_));
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* Esdc1a<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void Esdc1a<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;
        auto index     = [](Esdc1aInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::efd, [this, index]
                      { return y_.getData()[index(Esdc1aInternalVariables::EFD)]; });
        monitor_->set(Variable::vc, [this, index]
                      { return y_.getData()[index(Esdc1aInternalVariables::VC)]; });
        monitor_->set(Variable::vr, [this, index]
                      { return y_.getData()[index(Esdc1aInternalVariables::VR)]; });
        monitor_->set(Variable::vf, [this, index]
                      { return y_.getData()[index(Esdc1aInternalVariables::VF)]; });
        monitor_->set(Variable::se, [this, index]
                      { return y_.getData()[index(Esdc1aInternalVariables::SE)]; });
        monitor_->set(Variable::vfe, [this, index]
                      { return y_.getData()[index(Esdc1aInternalVariables::VFE)]; });
      }

      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::allocate()
      {
        size_     = static_cast<IdxT>(Esdc1aInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(2, ScalarT{0});

        auto signal_size = static_cast<size_t>(Esdc1aExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<Esdc1aInternalVariables::EFD>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<Esdc1aInternalVariables::EFD>()->set(
              &y[static_cast<size_t>(Esdc1aInternalVariables::EFD)],
              &(this->getVariableIndex(static_cast<IdxT>(Esdc1aInternalVariables::EFD))));
        }

        allocated_ = true;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Esdc1a: " << message << '\n';
            ret += 1;
          }
        };

        if (bus_ == nullptr)
        {
          Log::error() << "Esdc1a: bus pointer is null\n";
          ret += 1;
        }

        check(Ka_ > ZERO<RealT>, "Ka must be positive");
        check(Ta_ > ZERO<RealT>, "Ta must be positive");
        check(Tc_ >= ZERO<RealT>, "Tc must be non-negative");
        check(Te_ > ZERO<RealT>, "Te must be positive");
        check(Vrmin_ <= Vrmax_, "Vrmin must be less than or equal to Vrmax");
        check(spdmlt_ == ZERO<RealT> || spdmlt_ == ONE<RealT>, "Spdmlt must be 0 or 1");
        check(exclim_ == ZERO<RealT> || exclim_ == ONE<RealT>, "exclim must be 0 or 1");
        check(UEL_ >= static_cast<IdxT>(0) && UEL_ <= static_cast<IdxT>(3),
              "UEL must be 0, 1, 2, or 3");

        if (!(Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>) )
        {
          check(E1_ > ZERO<RealT>, "E1 must be positive when saturation is enabled");
          check(E2_ > ZERO<RealT>, "E2 must be positive when saturation is enabled");
          check(Se1_ > ZERO<RealT>, "Se1 must be positive when saturation is enabled");
          check(Se2_ > ZERO<RealT>, "Se2 must be positive when saturation is enabled");
          check(E1_ != E2_, "E1 and E2 must differ when saturation is enabled");
          check(Se1_ != Se2_, "Se1 and Se2 must differ when saturation is enabled");
        }

        if (!signals_.template isAssigned<Esdc1aInternalVariables::EFD>())
        {
          Log::error() << "Esdc1a: required EFD signal is not assigned\n";
          ret += 1;
        }

        if (spdmlt_ == ONE<RealT>
            && !signals_.template isAttached<Esdc1aExternalVariables::OMEGA>())
        {
          Log::error() << "Esdc1a: speed signal is required when Spdmlt is enabled\n";
          ret += 1;
        }

        auto check_attached_signal = [&](bool attached, bool linked, const char* name)
        {
          if (attached && !linked)
          {
            Log::error() << "Esdc1a: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_attached_signal(
            signals_.template isAttached<Esdc1aExternalVariables::OMEGA>(),
            signals_.template isAttached<Esdc1aExternalVariables::OMEGA>()
                && signals_.template isLinked<Esdc1aExternalVariables::OMEGA>(),
            "speed");
        check_attached_signal(
            signals_.template isAttached<Esdc1aExternalVariables::VREF>(),
            signals_.template isAttached<Esdc1aExternalVariables::VREF>()
                && signals_.template isLinked<Esdc1aExternalVariables::VREF>(),
            "VREF");
        check_attached_signal(
            signals_.template isAttached<Esdc1aExternalVariables::VS>(),
            signals_.template isAttached<Esdc1aExternalVariables::VS>()
                && signals_.template isLinked<Esdc1aExternalVariables::VS>(),
            "VS");
        check_attached_signal(
            signals_.template isAttached<Esdc1aExternalVariables::VUEL>(),
            signals_.template isAttached<Esdc1aExternalVariables::VUEL>()
                && signals_.template isLinked<Esdc1aExternalVariables::VUEL>(),
            "VUEL");

        return ret;
      }

      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::initialize()
      {
        if (verify() > 0)
        {
          Log::error() << "Esdc1a: cannot initialize with invalid configuration\n";
          return 1;
        }

        const auto EFDP = static_cast<size_t>(Esdc1aInternalVariables::EFDP);
        const auto VC   = static_cast<size_t>(Esdc1aInternalVariables::VC);
        const auto VR   = static_cast<size_t>(Esdc1aInternalVariables::VR);
        const auto VF   = static_cast<size_t>(Esdc1aInternalVariables::VF);
        const auto XLL  = static_cast<size_t>(Esdc1aInternalVariables::XLL);
        const auto EV   = static_cast<size_t>(Esdc1aInternalVariables::EV);
        const auto VLL  = static_cast<size_t>(Esdc1aInternalVariables::VLL);
        const auto VHV  = static_cast<size_t>(Esdc1aInternalVariables::VHV);
        const auto SE   = static_cast<size_t>(Esdc1aInternalVariables::SE);
        const auto VFE  = static_cast<size_t>(Esdc1aInternalVariables::VFE);
        const auto EFD  = static_cast<size_t>(Esdc1aInternalVariables::EFD);

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<Esdc1aExternalVariables::OMEGA>())
        {
          omega0 = signals_.template readExternalVariable<Esdc1aExternalVariables::OMEGA>();
        }

        ScalarT vs0{ZERO<RealT>};
        if (signals_.template isAttached<Esdc1aExternalVariables::VS>())
        {
          vs0 = signals_.template readExternalVariable<Esdc1aExternalVariables::VS>();
        }

        ScalarT vuel0{ZERO<RealT>};
        if (signals_.template isAttached<Esdc1aExternalVariables::VUEL>())
        {
          vuel0 = signals_.template readExternalVariable<Esdc1aExternalVariables::VUEL>();
        }

        const ScalarT d0 = ONE<RealT> + spdmlt_ * omega0;
        if (d0 == ZERO<RealT>)
        {
          Log::error() << "Esdc1a: speed multiplier denominator is zero at initialization\n";
          return 1;
        }

        const ScalarT Ec0 = std::sqrt(bus_->Vr() * bus_->Vr() + bus_->Vi() * bus_->Vi());

        const ScalarT efd0  = y[EFD];
        const ScalarT efdp0 = efd0 / d0;
        const ScalarT se0   = SB_ * Math::qramp(efdp0 - SA_);
        const ScalarT vfe0  = slim_off_ * (Ke_ + se0) * efdp0
                             + slim_ * Math::ramp((Ke_ + se0) * efdp0);
        const ScalarT vr0          = vfe0;
        const ScalarT vhv0         = vr0 / Ka_;
        auto          inverse_ramp = [](RealT y)
        {
          const RealT scaled_y = Math::MU<RealT> * y;
          if (scaled_y > static_cast<RealT>(50.0))
          {
            return y;
          }
          return std::log(std::expm1(scaled_y)) / Math::MU<RealT>;
        };

        ScalarT gate_input0 = vhv0;
        if (sUEL_ == ZERO<RealT>)
        {
          const RealT ramp_target = static_cast<RealT>(vhv0 - vuel0);
          if (ramp_target <= ZERO<RealT>)
          {
            Log::error() << "Esdc1a: smooth high-value gate is active at initialization\n";
            return 1;
          }
          gate_input0 = vuel0 + inverse_ramp(ramp_target);
        }

        const ScalarT vc0  = Ec0;
        const ScalarT vf0  = ScalarT{ZERO<RealT>};
        const ScalarT ev0  = gate_input0;
        const ScalarT xll0 = gate_input0;
        const ScalarT vll0 = gate_input0;

        if (vr0 < Vrmin_ || vr0 > Vrmax_)
        {
          Log::error() << "Esdc1a: initialized VR is outside limits\n";
          return 1;
        }

        vref_ = ev0 + vc0 + vf0 - vs0 - sUEL_ * vuel0;
        if (signals_.template isAttached<Esdc1aExternalVariables::VREF>())
        {
          signals_.template writeExternalVariable<Esdc1aExternalVariables::VREF>(vref_);
        }

        y[EFDP] = efdp0;
        y[VC]   = vc0;
        y[VR]   = vr0;
        y[VF]   = vf0;
        y[XLL]  = xll0;
        y[EV]   = ev0;
        y[VLL]  = vll0;
        y[VHV]  = vhv0;
        y[SE]   = se0;
        y[VFE]  = vfe0;
        y[EFD]  = efd0;

        for (IdxT i = 0; i < size_; ++i)
        {
          yp[i] = ZERO<RealT>;
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(Esdc1aInternalVariables::EFDP)] = true;
        tag_[static_cast<size_t>(Esdc1aInternalVariables::VC)]   = true;
        tag_[static_cast<size_t>(Esdc1aInternalVariables::VR)]   = true;
        tag_[static_cast<size_t>(Esdc1aInternalVariables::VF)]   = true;
        tag_[static_cast<size_t>(Esdc1aInternalVariables::XLL)]  = true;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Esdc1a<ScalarT, IdxT>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto EFDP = static_cast<size_t>(Esdc1aInternalVariables::EFDP);
        const auto VC   = static_cast<size_t>(Esdc1aInternalVariables::VC);
        const auto VR   = static_cast<size_t>(Esdc1aInternalVariables::VR);
        const auto VF   = static_cast<size_t>(Esdc1aInternalVariables::VF);
        const auto XLL  = static_cast<size_t>(Esdc1aInternalVariables::XLL);
        const auto EV   = static_cast<size_t>(Esdc1aInternalVariables::EV);
        const auto VLL  = static_cast<size_t>(Esdc1aInternalVariables::VLL);
        const auto VHV  = static_cast<size_t>(Esdc1aInternalVariables::VHV);
        const auto SE   = static_cast<size_t>(Esdc1aInternalVariables::SE);
        const auto VFE  = static_cast<size_t>(Esdc1aInternalVariables::VFE);
        const auto EFD  = static_cast<size_t>(Esdc1aInternalVariables::EFD);

        const auto OMEGA = static_cast<size_t>(Esdc1aExternalVariables::OMEGA);
        const auto VREF  = static_cast<size_t>(Esdc1aExternalVariables::VREF);
        const auto VS    = static_cast<size_t>(Esdc1aExternalVariables::VS);
        const auto VUEL  = static_cast<size_t>(Esdc1aExternalVariables::VUEL);

        const ScalarT efdp = y[EFDP];
        const ScalarT vc   = y[VC];
        const ScalarT vr   = y[VR];
        const ScalarT vf   = y[VF];
        const ScalarT xll  = y[XLL];
        const ScalarT ev   = y[EV];
        const ScalarT vll  = y[VLL];
        const ScalarT vhv  = y[VHV];
        const ScalarT se   = y[SE];
        const ScalarT vfe  = y[VFE];
        const ScalarT efd  = y[EFD];

        const ScalarT omega = ws[OMEGA];
        const ScalarT vref  = ws[VREF];
        const ScalarT vs    = ws[VS];
        const ScalarT vuel  = ws[VUEL];

        const ScalarT Ec        = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]);
        const ScalarT ev_target = vref + vs + sUEL_ * vuel - vc - vf;

        f[EFDP] = -yp[EFDP] + (vr - vfe) / Te_;
        f[VC]   = -yp[VC] + (Ec - vc) / Tr_;
        f[VR]   = -yp[VR] + Math::antiwindup(vr, -vr + Ka_ * vhv, Vrmin_, Vrmax_) / Ta_;
        f[VF]   = -yp[VF] + (-vf + Kf_ * (vr - vfe) / Te_) / Tf1_;
        f[XLL]  = -yp[XLL] + (ev - xll) / Tb_;
        f[EV]   = -ev + ev_target;
        f[VLL]  = -vll + xll + (Tc_ / Tb_) * (ev - xll);
        f[VHV]  = -vhv + sUEL_ * vll + sUELoff_ * Math::max(vll, vuel);
        f[SE]   = -se + SB_ * Math::qramp(efdp - SA_);
        f[VFE]  = -vfe + slim_off_ * (Ke_ + se) * efdp + slim_ * Math::ramp((Ke_ + se) * efdp);
        f[EFD]  = -efd + (ONE<RealT> + spdmlt_ * omega) * efdp;

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::evaluateResidual()
      {
        const auto OMEGA = static_cast<size_t>(Esdc1aExternalVariables::OMEGA);
        const auto VREF  = static_cast<size_t>(Esdc1aExternalVariables::VREF);
        const auto VS    = static_cast<size_t>(Esdc1aExternalVariables::VS);
        const auto VUEL  = static_cast<size_t>(Esdc1aExternalVariables::VUEL);

        std::fill(ws_.begin(), ws_.end(), ScalarT{ZERO<RealT>});
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);
        ws_[VREF] = vref_;

        if (signals_.template isAttached<Esdc1aExternalVariables::OMEGA>())
        {
          ws_[OMEGA] = signals_.template readExternalVariable<Esdc1aExternalVariables::OMEGA>();
          ws_indices_[OMEGA] =
              signals_.template readExternalVariableIndex<Esdc1aExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<Esdc1aExternalVariables::VREF>())
        {
          ws_[VREF] = signals_.template readExternalVariable<Esdc1aExternalVariables::VREF>();
          ws_indices_[VREF] =
              signals_.template readExternalVariableIndex<Esdc1aExternalVariables::VREF>();
        }
        if (signals_.template isAttached<Esdc1aExternalVariables::VS>())
        {
          ws_[VS] = signals_.template readExternalVariable<Esdc1aExternalVariables::VS>();
          ws_indices_[VS] =
              signals_.template readExternalVariableIndex<Esdc1aExternalVariables::VS>();
        }
        if (signals_.template isAttached<Esdc1aExternalVariables::VUEL>())
        {
          ws_[VUEL] = signals_.template readExternalVariable<Esdc1aExternalVariables::VUEL>();
          ws_indices_[VUEL] =
              signals_.template readExternalVariableIndex<Esdc1aExternalVariables::VUEL>();
        }

        wb_[0] = bus_->Vr();
        wb_[1] = bus_->Vi();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);

        f_.setDataUpdated();
        return 0;
      }
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit

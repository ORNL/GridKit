/**
 * @file GastPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the GASTPTI governor model.
 */

#pragma once

#include <algorithm>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::GastPti()
      {
        size_ = static_cast<IdxT>(GastPtiInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::GastPti(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(GastPtiInternalVariables::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::~GastPti()
      {
      }

      template <typename scalar_type, typename index_type>
      scalar_type GastPti<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * va_system_base_ / va_component_base_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type GastPti<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<scalar_type>(ONE<RealT>));
      }

      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        if (data.parameters.contains(Params::R))
        {
          R_ = std::get<RealT>(data.parameters.at(Params::R));
        }
        if (data.parameters.contains(Params::T1))
        {
          T1_ = std::get<RealT>(data.parameters.at(Params::T1));
        }
        if (data.parameters.contains(Params::T2))
        {
          T2_ = std::get<RealT>(data.parameters.at(Params::T2));
        }
        if (data.parameters.contains(Params::T3))
        {
          T3_ = std::get<RealT>(data.parameters.at(Params::T3));
        }
        if (data.parameters.contains(Params::At))
        {
          At_ = std::get<RealT>(data.parameters.at(Params::At));
        }
        if (data.parameters.contains(Params::Kt))
        {
          Kt_ = std::get<RealT>(data.parameters.at(Params::Kt));
        }
        if (data.parameters.contains(Params::Vmax))
        {
          Vmax_ = std::get<RealT>(data.parameters.at(Params::Vmax));
        }
        if (data.parameters.contains(Params::Vmin))
        {
          Vmin_ = std::get<RealT>(data.parameters.at(Params::Vmin));
        }
        if (data.parameters.contains(Params::Dturb))
        {
          Dturb_ = std::get<RealT>(data.parameters.at(Params::Dturb));
        }
        if (data.parameters.contains(Params::Trate))
        {
          Trate_ = std::get<RealT>(data.parameters.at(Params::Trate));
        }
        if (data.parameters.contains(Params::mode))
        {
          mode_ = static_cast<ResponseMode>(
              std::get<IdxT>(data.parameters.at(Params::mode)));
        }

        if (mode_ == ResponseMode::DownOnly)
        {
          Log::warning() << "GastPti: mode 2 (Down Only) is not yet supported; "
                            "using mode 0 (Normal)\n";
        }

        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::setDerivedParameters()
      {
        va_component_base_ = Trate_ * static_cast<RealT>(1.0e6);

        T1_ = std::max(T1_, TIME_CONSTANT_MINIMUM);
        T2_ = std::max(T2_, TIME_CONSTANT_MINIMUM);
        T3_ = std::max(T3_, TIME_CONSTANT_MINIMUM);

        sfix_ = ONE<RealT>;
        if (mode_ == ResponseMode::Fixed)
        {
          sfix_ = ZERO<RealT>;
        }
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* GastPti<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        auto index     = [](GastPtiInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::pmech, [this, index]
                      { return y_.getData()[index(GastPtiInternalVariables::PMECH)]; });
        monitor_->set(Variable::xvalve, [this, index]
                      { return y_.getData()[index(GastPtiInternalVariables::XVALVE)]; });
        monitor_->set(Variable::xflow, [this, index]
                      { return y_.getData()[index(GastPtiInternalVariables::XFLOW)]; });
        monitor_->set(Variable::xtemp, [this, index]
                      { return y_.getData()[index(GastPtiInternalVariables::XTEMP)]; });
        monitor_->set(Variable::vload, [this, index]
                      { return y_.getData()[index(GastPtiInternalVariables::VLOAD)]; });
        monitor_->set(Variable::vtemp, [this, index]
                      { return y_.getData()[index(GastPtiInternalVariables::VTEMP)]; });
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.clear();

        auto signal_size = static_cast<size_t>(GastPtiExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<GastPtiInternalVariables::PMECH>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<GastPtiInternalVariables::PMECH>()->set(
              &y[static_cast<size_t>(GastPtiInternalVariables::PMECH)],
              &(this->getVariableIndex(static_cast<IdxT>(GastPtiInternalVariables::PMECH))));
        }

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::verify() const
      {
        int ret = 0;

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "GastPti: " << message << '\n';
            ret += 1;
          }
        };

        check(R_ > ZERO<RealT>, "R must be positive");
        check(At_ >= ZERO<RealT>, "At must be non-negative");
        check(Kt_ >= ZERO<RealT>, "Kt must be non-negative");
        const bool valid_mode = mode_ == ResponseMode::Normal
                                || mode_ == ResponseMode::Fixed
                                || mode_ == ResponseMode::DownOnly;
        check(valid_mode, "mode must be 0 (Normal), 1 (Fixed), or 2 (Down Only)");

        if (mode_ == ResponseMode::Fixed)
        {
          check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax in Fixed mode");
        }
        else
        {
          check(Vmin_ < Vmax_, "Vmin must be less than Vmax in Normal mode");
        }

        if (mode_ != ResponseMode::Fixed && Vmin_ < Vmax_)
        {
          constexpr RealT minimum_gate{0.99};

          const RealT midpoint = HALF<RealT> * (Vmin_ + Vmax_);
          const RealT maximum_gate =
              Math::indicator(midpoint, ZERO<RealT>, Vmin_, Vmax_);

          if (maximum_gate < minimum_gate)
          {
            Log::warning()
                << "GastPti: narrow valve limits reduce the maximum anti-windup gate to "
                << maximum_gate << "; Normal-mode dynamics may be attenuated\n";
          }
        }

        check(Dturb_ >= ZERO<RealT>, "Dturb must be non-negative");
        check(Trate_ > ZERO<RealT>, "Trate must be positive");

        if (signals_.template isAttached<GastPtiExternalVariables::OMEGA>()
            && !signals_.template isLinked<GastPtiExternalVariables::OMEGA>())
        {
          Log::error() << "GastPti: omega signal attached with no linked source\n";
          ret += 1;
        }

        if (signals_.template isAttached<GastPtiExternalVariables::PREF>()
            && !signals_.template isLinked<GastPtiExternalVariables::PREF>())
        {
          Log::error() << "GastPti: pref signal attached with no linked source\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::initialize()
      {
        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        const auto XVALVE = static_cast<size_t>(GastPtiInternalVariables::XVALVE);
        const auto XFLOW  = static_cast<size_t>(GastPtiInternalVariables::XFLOW);
        const auto XTEMP  = static_cast<size_t>(GastPtiInternalVariables::XTEMP);
        const auto VLOAD  = static_cast<size_t>(GastPtiInternalVariables::VLOAD);
        const auto VTEMP  = static_cast<size_t>(GastPtiInternalVariables::VTEMP);
        const auto VLV    = static_cast<size_t>(GastPtiInternalVariables::VLV);
        const auto PMECH  = static_cast<size_t>(GastPtiInternalVariables::PMECH);

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<GastPtiExternalVariables::OMEGA>())
        {
          omega0 = signals_.template readExternalVariable<GastPtiExternalVariables::OMEGA>();
        }

        // Machine initialization seeds the assigned PMECH signal before
        // governor initialization reads it.
        const ScalarT pmech0 = toComponentBase(y[PMECH]);
        const ScalarT xflow0 = pmech0 + Dturb_ * omega0;
        const ScalarT vtemp0 = At_ + Kt_ * (At_ - xflow0);

        if (mode_ != ResponseMode::Fixed
            && (static_cast<RealT>(xflow0) < Vmin_ || static_cast<RealT>(xflow0) > Vmax_))
        {
          Log::warning() << "GastPti: initialization of governor from power flow exceeded rating\n";
        }

        const RealT margin = static_cast<RealT>(vtemp0 - xflow0);
        if (margin <= ZERO<RealT>)
        {
          Log::error() << "GastPti: initial temperature-gate margin must be positive\n";
          return 1;
        }

        // The temperature gate is the smooth minimum, so VLOAD must sit where
        // min(VLOAD, VTEMP) reproduces the initial fuel flow exactly.
        const ScalarT vload0 = vtemp0 - Math::iramp(margin);

        y[XVALVE] = xflow0;
        y[XFLOW]  = xflow0;
        y[XTEMP]  = xflow0;
        y[VLOAD]  = vload0;
        y[VTEMP]  = vtemp0;
        y[VLV]    = xflow0;
        y[PMECH]  = toSystemBase(pmech0);

        pref_set_ = toSystemBase(vload0 + omega0 / R_);
        if (signals_.template isAttached<GastPtiExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<GastPtiExternalVariables::PREF>(pref_set_);
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
      int GastPti<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(GastPtiInternalVariables::XVALVE)] = true;
        tag_[static_cast<size_t>(GastPtiInternalVariables::XFLOW)]  = true;
        tag_[static_cast<size_t>(GastPtiInternalVariables::XTEMP)]  = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int GastPti<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        const auto XVALVE = static_cast<size_t>(GastPtiInternalVariables::XVALVE);
        const auto XFLOW  = static_cast<size_t>(GastPtiInternalVariables::XFLOW);
        const auto XTEMP  = static_cast<size_t>(GastPtiInternalVariables::XTEMP);
        const auto VLOAD  = static_cast<size_t>(GastPtiInternalVariables::VLOAD);
        const auto VTEMP  = static_cast<size_t>(GastPtiInternalVariables::VTEMP);
        const auto VLV    = static_cast<size_t>(GastPtiInternalVariables::VLV);
        const auto PMECH  = static_cast<size_t>(GastPtiInternalVariables::PMECH);

        const auto OMEGA = static_cast<size_t>(GastPtiExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(GastPtiExternalVariables::PREF);

        const ScalarT xvalve = y[XVALVE];
        const ScalarT xflow  = y[XFLOW];
        const ScalarT xtemp  = y[XTEMP];
        const ScalarT vload  = y[VLOAD];
        const ScalarT vtemp  = y[VTEMP];
        const ScalarT vlv    = y[VLV];
        const ScalarT pmech  = y[PMECH];

        const ScalarT omega = ws[OMEGA];
        const ScalarT pref  = toComponentBase(ws[PREF]);

        const ScalarT fvalve       = vlv - xvalve;
        const ScalarT fflow        = -xflow + xvalve;
        const ScalarT ftemp        = -xtemp + xflow;
        const ScalarT valve_target = Math::antiwindup(xvalve, fvalve, Vmin_, Vmax_);

        f[XVALVE] = -yp[XVALVE] + sfix_ * valve_target / T1_;
        f[XFLOW]  = -yp[XFLOW] + sfix_ * fflow / T2_;
        f[XTEMP]  = -yp[XTEMP] + sfix_ * ftemp / T3_;
        f[VLOAD]  = -omega + R_ * (pref - vload);
        f[VTEMP]  = -vtemp + At_ + Kt_ * (At_ - xtemp);
        f[VLV]    = -vlv + Math::min(vload, vtemp);
        f[PMECH]  = -toComponentBase(pmech) + xflow - Dturb_ * omega;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::evaluateResidual()
      {
        const auto OMEGA = static_cast<size_t>(GastPtiExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(GastPtiExternalVariables::PREF);

        ws_[OMEGA]         = ZERO<RealT>;
        ws_[PREF]          = pref_set_;
        ws_indices_[OMEGA] = INVALID_INDEX<IdxT>;
        ws_indices_[PREF]  = INVALID_INDEX<IdxT>;

        if (signals_.template isAttached<GastPtiExternalVariables::OMEGA>())
        {
          ws_[OMEGA] = signals_.template readExternalVariable<GastPtiExternalVariables::OMEGA>();
          ws_indices_[OMEGA] =
              signals_.template readExternalVariableIndex<GastPtiExternalVariables::OMEGA>();
        }

        if (signals_.template isAttached<GastPtiExternalVariables::PREF>())
        {
          ws_[PREF] = signals_.template readExternalVariable<GastPtiExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<GastPtiExternalVariables::PREF>();
        }

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);
        f_.setDataUpdated();
        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

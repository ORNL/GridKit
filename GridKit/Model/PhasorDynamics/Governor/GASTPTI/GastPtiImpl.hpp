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

      template <class ScalarT, typename IdxT>
      GastPti<ScalarT, IdxT>::GastPti()
      {
        size_ = static_cast<IdxT>(GastPtiInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      GastPti<ScalarT, IdxT>::GastPti(const model_data_type& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(GastPtiInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      GastPti<ScalarT, IdxT>::~GastPti()
      {
      }

      template <class ScalarT, typename IdxT>
      void GastPti<ScalarT, IdxT>::initModelParams(const model_data_type& data)
      {
        using Params = typename model_data_type::Parameters;

        parameter_error_count_ = 0;

        auto load_required_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "GastPti: missing required parameter '" << name << "'\n";
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
            Log::error() << "GastPti: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto load_optional_real = [&](auto key, RealT& target, const char* name)
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
            Log::error() << "GastPti: optional parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        load_required_real(Params::R, R_, "R");
        load_required_real(Params::T1, T1_, "T1");
        load_required_real(Params::T2, T2_, "T2");
        load_required_real(Params::T3, T3_, "T3");
        load_required_real(Params::At, At_, "At");
        load_required_real(Params::Kt, Kt_, "Kt");
        load_required_real(Params::Vmax, Vmax_, "Vmax");
        load_required_real(Params::Vmin, Vmin_, "Vmin");
        load_required_real(Params::Dturb, Dturb_, "Dturb");
        load_optional_real(Params::Trate, Trate_, "Trate");
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* GastPti<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void GastPti<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;
        auto index     = [](GastPtiInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::pmech, [this, index]
                      { return y_[index(GastPtiInternalVariables::PMECH)]; });
        monitor_->set(Variable::fuelvalve, [this, index]
                      { return y_[index(GastPtiInternalVariables::XVALVE)]; });
        monitor_->set(Variable::fuelflow, [this, index]
                      { return y_[index(GastPtiInternalVariables::XFLOW)]; });
        monitor_->set(Variable::exhausttemp, [this, index]
                      { return y_[index(GastPtiInternalVariables::XTEMP)]; });
        monitor_->set(Variable::vload, [this, index]
                      { return y_[index(GastPtiInternalVariables::VLOAD)]; });
        monitor_->set(Variable::vtemp, [this, index]
                      { return y_[index(GastPtiInternalVariables::VTEMP)]; });
        monitor_->set(Variable::vlv, [this, index]
                      { return y_[index(GastPtiInternalVariables::VLV)]; });
      }

      template <class ScalarT, typename IdxT>
      int GastPti<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int GastPti<ScalarT, IdxT>::allocate()
      {
        size_     = static_cast<IdxT>(GastPtiInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        f_.assign(size, ScalarT{0});
        y_.assign(size, ScalarT{0});
        yp_.assign(size, ScalarT{0});
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
          signals_.template getSignalNode<GastPtiInternalVariables::PMECH>()->set(
              &y_[static_cast<size_t>(GastPtiInternalVariables::PMECH)],
              &(this->getVariableIndex(static_cast<IdxT>(GastPtiInternalVariables::PMECH))));
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int GastPti<ScalarT, IdxT>::verify() const
      {
        int ret = parameter_error_count_;

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "GastPti: " << message << '\n';
            ret += 1;
          }
        };

        check(R_ > ZERO<RealT>, "R must be positive");
        check(T1_ >= ZERO<RealT>, "T1 must be non-negative");
        check(T2_ >= ZERO<RealT>, "T2 must be non-negative");
        check(T3_ >= ZERO<RealT>, "T3 must be non-negative");
        check(At_ >= ZERO<RealT>, "At must be non-negative");
        check(Kt_ >= ZERO<RealT>, "Kt must be non-negative");
        check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax");
        check(Dturb_ >= ZERO<RealT>, "Dturb must be non-negative");
        check(Trate_ >= ZERO<RealT>, "Trate must be non-negative");

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

      template <class ScalarT, typename IdxT>
      int GastPti<ScalarT, IdxT>::initialize()
      {
        if (parameter_error_count_ > 0 || verify() > 0)
        {
          Log::error() << "GastPti: cannot initialize with invalid configuration\n";
          return 1;
        }

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
        const ScalarT pmech0 = y_[PMECH];
        const ScalarT xflow0 = pmech0 + Dturb_ * omega0;

        y_[XFLOW]  = xflow0;
        y_[XVALVE] = xflow0;
        y_[XTEMP]  = xflow0;
        y_[VTEMP]  = At_ + Kt_ * (At_ - y_[XTEMP]);
        y_[VLV]    = y_[XVALVE];
        y_[VLOAD]  = y_[XVALVE];
        y_[PMECH]  = pmech0;

        pref_set_ = y_[VLOAD] + omega0 / R_;
        if (signals_.template isAttached<GastPtiExternalVariables::PREF>())
        {
          pref_set_ = signals_.template readExternalVariable<GastPtiExternalVariables::PREF>();
          y_[VLOAD] = pref_set_ - omega0 / R_;
        }

        std::fill(yp_.begin(), yp_.end(), ZERO<RealT>);
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int GastPti<ScalarT, IdxT>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(GastPtiInternalVariables::XVALVE)] = (T1_ > ZERO<RealT>);
        tag_[static_cast<size_t>(GastPtiInternalVariables::XFLOW)]  = (T2_ > ZERO<RealT>);
        tag_[static_cast<size_t>(GastPtiInternalVariables::XTEMP)]  = (T3_ > ZERO<RealT>);
        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int GastPti<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT*                  y,
          ScalarT*                  yp,
          [[maybe_unused]] ScalarT* wb,
          ScalarT*                  ws,
          ScalarT*                  f)
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
        const ScalarT pref  = ws[PREF];

        const ScalarT fvalve = vlv - xvalve;

        f[XVALVE] = -T1_ * yp[XVALVE] + Math::antiwindup(xvalve, fvalve, Vmin_, Vmax_);
        f[XFLOW]  = -T2_ * yp[XFLOW] - xflow + xvalve;
        f[XTEMP]  = -T3_ * yp[XTEMP] - xtemp + xflow;
        f[VLOAD]  = -vload + pref - omega / R_;
        f[VTEMP]  = -vtemp + At_ + Kt_ * (At_ - xtemp);
        f[VLV]    = -vlv + Math::min(vload, vtemp);
        f[PMECH]  = -pmech + xflow - Dturb_ * omega;

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int GastPti<ScalarT, IdxT>::evaluateResidual()
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

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit

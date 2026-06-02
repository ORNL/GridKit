/**
 * @file ForcedOscillationImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the forced oscillation signal operator.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <variant>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalOperator/ForcedOscillation/ForcedOscillation.hpp>
#include <GridKit/Model/PhasorDynamics/SignalOperator/ForcedOscillation/ForcedOscillationData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalOperator
    {
      using Log = ::GridKit::Utilities::Logger;

      template <class ScalarT, typename IdxT>
      ForcedOscillation<ScalarT, IdxT>::ForcedOscillation()
      {
        size_  = 1;
        time_  = ZERO<RealT>;
        alpha_ = ZERO<RealT>;
      }

      template <class ScalarT, typename IdxT>
      ForcedOscillation<ScalarT, IdxT>::ForcedOscillation(const model_data_type& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_  = 1;
        time_  = ZERO<RealT>;
        alpha_ = ZERO<RealT>;
      }

      template <class ScalarT, typename IdxT>
      ForcedOscillation<ScalarT, IdxT>::~ForcedOscillation()
      {
      }

      template <class ScalarT, typename IdxT>
      typename ForcedOscillation<ScalarT, IdxT>::RealT
      ForcedOscillation<ScalarT, IdxT>::readRealParameter(
          const model_data_type&               data,
          typename model_data_type::Parameters key,
          RealT                                fallback)
      {
        if (!data.parameters.contains(key))
        {
          return fallback;
        }

        const auto& value = data.parameters.at(key);
        if (const auto* real_value = std::get_if<RealT>(&value))
        {
          return *real_value;
        }
        if (const auto* index_value = std::get_if<IdxT>(&value))
        {
          return static_cast<RealT>(*index_value);
        }

        Log::error() << "ForcedOscillation: parameter "
                     << static_cast<int>(key)
                     << " must be numeric\n";
        parameter_error_count_ += 1;
        return fallback;
      }

      template <class ScalarT, typename IdxT>
      void ForcedOscillation<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
      {
        using Parameters = typename model_data_type::Parameters;

        A_       = readRealParameter(data, Parameters::A, A_);
        f_param_ = readRealParameter(data, Parameters::f, f_param_);
        Kf_      = readRealParameter(data, Parameters::Kf, Kf_);
        Phi_     = readRealParameter(data, Parameters::Phi, Phi_);
        Bias_    = readRealParameter(data, Parameters::Bias, Bias_);
        Kin_     = readRealParameter(data, Parameters::Kin, Kin_);
        u0_      = readRealParameter(data, Parameters::u0, u0_);
        Ton_     = readRealParameter(data, Parameters::Ton, Ton_);
        Toff_    = readRealParameter(data, Parameters::Toff, Toff_);
        Tr_      = readRealParameter(data, Parameters::Tr, Tr_);
        Tf_      = readRealParameter(data, Parameters::Tf, Tf_);
        Kd_      = readRealParameter(data, Parameters::Kd, Kd_);
        Lmin_    = readRealParameter(data, Parameters::Lmin, Lmin_);
        Lmax_    = readRealParameter(data, Parameters::Lmax, Lmax_);
      }

      template <class ScalarT, typename IdxT>
      int ForcedOscillation<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int ForcedOscillation<ScalarT, IdxT>::allocate()
      {
        auto size = static_cast<size_t>(size_);

        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(1, ScalarT{0.0});
        h_.clear();
        ws_.assign(1, ScalarT{0.0});
        ws_indices_.assign(1, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<ForcedOscillationInternalVariables::OUT>())
        {
          signals_.template getSignalNode<ForcedOscillationInternalVariables::OUT>()->set(
              &y_[0], &(this->getVariableIndex(0)));
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int ForcedOscillation<ScalarT, IdxT>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "ForcedOscillation: " << message << '\n';
            ret += 1;
          }
        };

        if (!signals_.template isAssigned<ForcedOscillationInternalVariables::OUT>())
        {
          Log::error() << "ForcedOscillation: required output signal is not assigned\n";
          ret += 1;
        }

        if (signals_.template isAttached<ForcedOscillationExternalVariables::IN>()
            && !signals_.template isLinked<ForcedOscillationExternalVariables::IN>())
        {
          Log::error() << "ForcedOscillation: input signal attached with no linked source\n";
          ret += 1;
        }

        check(A_ >= ZERO<RealT>, "A must be non-negative");
        check(f_param_ >= ZERO<RealT>, "f must be non-negative");
        check(Kf_ >= ZERO<RealT>, "Kf must be non-negative");
        check(Tr_ >= ZERO<RealT>, "Tr must be non-negative");
        check(Tf_ >= ZERO<RealT>, "Tf must be non-negative");
        check(Kd_ >= ZERO<RealT>, "Kd must be non-negative");
        check(Toff_ < ZERO<RealT> || Toff_ >= Ton_, "Toff must be negative or greater than or equal to Ton");
        check(Lmin_ <= Lmax_, "Lmin must be less than or equal to Lmax");

        return ret;
      }

      template <class ScalarT, typename IdxT>
      int ForcedOscillation<ScalarT, IdxT>::initialize()
      {
        refreshForcing(time_);
        refreshInitialOutput();

        yp_[0] = ZERO<RealT>;

        return 0;
      }

      template <class ScalarT, typename IdxT>
      void ForcedOscillation<ScalarT, IdxT>::updateTime(RealT t, RealT a)
      {
        time_  = t;
        alpha_ = a;

        refreshForcing(time_);

        in_ = u0_;
        if (signals_.template isAttached<ForcedOscillationExternalVariables::IN>()
            && signals_.template isLinked<ForcedOscillationExternalVariables::IN>())
        {
          in_ = signals_.template readExternalVariable<ForcedOscillationExternalVariables::IN>();
        }
        vraw_ = Kin_ * in_ + Bias_ + force_;
      }

      template <class ScalarT, typename IdxT>
      void ForcedOscillation<ScalarT, IdxT>::initializeInputFromOutput()
      {
        if (!signals_.template isAttached<ForcedOscillationExternalVariables::IN>()
            || !signals_.template isLinked<ForcedOscillationExternalVariables::IN>())
        {
          return;
        }

        if (Kin_ == ZERO<RealT>)
        {
          return;
        }

        refreshForcing(time_);

        const ScalarT output = y_[0];
        in_                  = (output - Bias_ - force_) / Kin_;
        vraw_                = Kin_ * in_ + Bias_ + force_;
        signals_.template writeExternalVariable<ForcedOscillationExternalVariables::IN>(in_);
      }

      template <class ScalarT, typename IdxT>
      int ForcedOscillation<ScalarT, IdxT>::tagDifferentiable()
      {
        tag_[0] = false;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      void ForcedOscillation<ScalarT, IdxT>::refreshForcing(RealT time)
      {
        static constexpr RealT pi     = static_cast<RealT>(3.141592653589793238462643383279502884L);
        static constexpr RealT two_pi = static_cast<RealT>(2.0) * pi;

        bool is_active = (time >= Ton_);
        if (Toff_ >= ZERO<RealT> && time >= Toff_)
        {
          is_active = false;
        }

        RealT envelope = ZERO<RealT>;
        if (is_active)
        {
          envelope = ONE<RealT>;

          if (Tr_ > ZERO<RealT> && time < Ton_ + Tr_)
          {
            const RealT x  = (time - Ton_) / Tr_;
            envelope      *= HALF<RealT> * (ONE<RealT> - std::cos(pi * x));
          }

          if (Toff_ >= ZERO<RealT> && Tf_ > ZERO<RealT> && time > Toff_ - Tf_)
          {
            const RealT x  = (Toff_ - time) / Tf_;
            envelope      *= HALF<RealT> * (ONE<RealT> - std::cos(pi * x));
          }
        }

        const RealT tau   = std::max(time - Ton_, ZERO<RealT>);
        const RealT phase = Phi_ + two_pi * (f_param_ * tau + HALF<RealT> * Kf_ * tau * tau);
        const RealT decay = std::exp(-Kd_ * tau);

        active_ = static_cast<RealT>(is_active ? ONE<RealT> : ZERO<RealT>);
        env_    = envelope;
        force_  = A_ * envelope * decay * std::sin(phase);
      }

      template <class ScalarT, typename IdxT>
      void ForcedOscillation<ScalarT, IdxT>::refreshInitialOutput()
      {
        in_ = u0_;
        if (signals_.template isAttached<ForcedOscillationExternalVariables::IN>()
            && signals_.template isLinked<ForcedOscillationExternalVariables::IN>())
        {
          in_ = signals_.template readExternalVariable<ForcedOscillationExternalVariables::IN>();
        }

        vraw_ = Kin_ * in_ + Bias_ + force_;
        y_[0] = Math::clamp(vraw_, Lmin_, Lmax_);
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int ForcedOscillation<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT*                  y,
          [[maybe_unused]] ScalarT* yp,
          ScalarT*                  wb,
          ScalarT*                  ws,
          ScalarT*                  f)
      {
        ScalarT raw = Kin_ * ws[0] + wb[0];
        ScalarT out = Math::clamp(raw, Lmin_, Lmax_);

        f[0] = y[0] - out;

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int ForcedOscillation<ScalarT, IdxT>::evaluateResidual()
      {
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;

        bool input_attached = signals_.template isAttached<ForcedOscillationExternalVariables::IN>();
        if (input_attached)
        {
          ws_[0] = signals_.template readExternalVariable<ForcedOscillationExternalVariables::IN>();
          ws_indices_[0] =
              signals_.template readExternalVariableIndex<ForcedOscillationExternalVariables::IN>();
        }

        refreshForcing(time_);

        wb_[0] = Bias_ + force_;
        if (!input_attached)
        {
          wb_[0] += Kin_ * u0_;
        }

        in_ = u0_;
        if (input_attached)
        {
          in_ = ws_[0];
        }

        vraw_ = Kin_ * in_ + Bias_ + force_;

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

        return 0;
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* ForcedOscillation<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void ForcedOscillation<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;

        monitor_->set(Variable::in, [this]
                      { return in_; });
        monitor_->set(Variable::env, [this]
                      { return env_; });
        monitor_->set(Variable::force, [this]
                      { return force_; });
        monitor_->set(Variable::vraw, [this]
                      { return vraw_; });
        monitor_->set(Variable::out, [this]
                      { return y_[0]; });
        monitor_->set(Variable::active, [this]
                      { return active_; });
      }

    } // namespace SignalOperator
  } // namespace PhasorDynamics
} // namespace GridKit

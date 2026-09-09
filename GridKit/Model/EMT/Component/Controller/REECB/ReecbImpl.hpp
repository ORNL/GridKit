#pragma once

#include <cmath>
#include <mutex>
#include <numbers>

#include <GridKit/Model/EMT/Component/Controller/REECB/Reecb.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::Reecb()
        : Reecb(ModelDataT{})
      {
      }

      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::Reecb(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        size_ = static_cast<IdxT>(ReecbInternalVariables::MAXIMUM);
        signals_.template assignSignal<ReecbInternalVariables::ICMDD>(&output_[0]);
        signals_.template assignSignal<ReecbInternalVariables::ICMDQ>(&output_[1]);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::~Reecb() = default;

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;
        S_              = parameter<RealT>(data, Parameter::S, S_);
        V_              = parameter<RealT>(data, Parameter::V, V_);
        PfFlag_         = parameter<bool>(data, Parameter::PfFlag, PfFlag_);
        VFlag_          = parameter<bool>(data, Parameter::VFlag, VFlag_);
        QFlag_          = parameter<bool>(data, Parameter::QFlag, QFlag_);
        Pqflag_         = parameter<bool>(data, Parameter::Pqflag, Pqflag_);
        Trv_            = parameter<RealT>(data, Parameter::Trv, Trv_);
        Tp_             = parameter<RealT>(data, Parameter::Tp, Tp_);
        Vref0_          = parameter<RealT>(data, Parameter::Vref0, Vref0_);
        Vdip_           = parameter<RealT>(data, Parameter::Vdip, Vdip_);
        Vup_            = parameter<RealT>(data, Parameter::Vup, Vup_);
        dbd1_           = parameter<RealT>(data, Parameter::dbd1, dbd1_);
        dbd2_           = parameter<RealT>(data, Parameter::dbd2, dbd2_);
        kqv_            = parameter<RealT>(data, Parameter::kqv, kqv_);
        Iql1_           = parameter<RealT>(data, Parameter::Iql1, Iql1_);
        Iqh1_           = parameter<RealT>(data, Parameter::Iqh1, Iqh1_);
        Qmax_           = parameter<RealT>(data, Parameter::Qmax, Qmax_);
        Qmin_           = parameter<RealT>(data, Parameter::Qmin, Qmin_);
        Kqp_            = parameter<RealT>(data, Parameter::Kqp, Kqp_);
        Kqi_            = parameter<RealT>(data, Parameter::Kqi, Kqi_);
        Vmax_           = parameter<RealT>(data, Parameter::Vmax, Vmax_);
        Vmin_           = parameter<RealT>(data, Parameter::Vmin, Vmin_);
        Kvp_            = parameter<RealT>(data, Parameter::Kvp, Kvp_);
        Kvi_            = parameter<RealT>(data, Parameter::Kvi, Kvi_);
        Tiq_            = parameter<RealT>(data, Parameter::Tiq, Tiq_);
        Tpord_          = parameter<RealT>(data, Parameter::Tpord, Tpord_);
        dPmax_          = parameter<RealT>(data, Parameter::dPmax, dPmax_);
        dPmin_          = parameter<RealT>(data, Parameter::dPmin, dPmin_);
        Pmax_           = parameter<RealT>(data, Parameter::Pmax, Pmax_);
        Pmin_           = parameter<RealT>(data, Parameter::Pmin, Pmin_);
        Imax_           = parameter<RealT>(data, Parameter::Imax, Imax_);
        Vref0_given_    = data.parameters.contains(Parameter::Vref0);
        for (const auto value : {Trv_, Tp_, Tiq_, Tpord_})
          if (!std::isfinite(value) || value < ZERO<RealT>)
            throw std::invalid_argument("Reecb: time constants must be finite and nonnegative");
        if (std::min({Trv_, Tp_, Tiq_, Tpord_}) < TIME_CONSTANT_MINIMUM)
        {
          static std::once_flag warning;
          std::call_once(warning, []
                         { Log::warning() << "Reecb: Trv, Tp, Tiq, and Tpord below 0.001 s are raised to that floor\n"; });
        }
        Trv_      = std::max(Trv_, TIME_CONSTANT_MINIMUM);
        Tp_       = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        Tiq_      = std::max(Tiq_, TIME_CONSTANT_MINIMUM);
        Tpord_    = std::max(Tpord_, TIME_CONSTANT_MINIMUM);
        pf_on_    = static_cast<RealT>(PfFlag_);
        pf_off_   = ONE<RealT> - pf_on_;
        q_on_     = static_cast<RealT>(QFlag_);
        q_off_    = ONE<RealT> - q_on_;
        q_pi_on_  = static_cast<RealT>(QFlag_ && VFlag_);
        v_ref_on_ = static_cast<RealT>(QFlag_ && !VFlag_);
        q_ref_on_ = ONE<RealT> - v_ref_on_;
        pq_on_    = static_cast<RealT>(Pqflag_);
        pq_off_   = ONE<RealT> - pq_on_;
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::attachInput(InputSignals inputs)
      {
        if (allocated_)
          throw std::logic_error("Reecb inputs cannot change after allocation");
        for (size_t p = 0; p < inputs.size(); ++p)
          if (inputs[p])
            signals_.attachSignal(static_cast<ReecbExternalVariables>(p), inputs[p]);
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        if (allocated_)
          throw std::logic_error("Reecb outputs cannot change after allocation");
        const auto n = static_cast<size_t>(output);
        if (n >= output_.size() || !signal || alias_[n])
          throw std::invalid_argument("Reecb: invalid or duplicate output assignment");
        signal->claimProducer();
        alias_[n] = signal;
      }

      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::SignalT& Reecb<scalar_type, index_type>::outputSignal(Outputs output)
      {
        return output_.at(static_cast<size_t>(output));
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
          this->allocateVectors(size_);
        tag_.resize(static_cast<size_t>(size_));
        variable_indices_.resize(static_cast<size_t>(size_));
        residual_indices_.resize(static_cast<size_t>(size_));
        this->assignGlobalIndices(0);
        this->allocateExternalVectors(static_cast<IdxT>(ReecbExternalVariables::MAXIMUM), 0);
        signals_.registerExternalVariableSignals(*this);
        signals_.bindInternalVariableSignals(*this);
        for (size_t n = 0; n < alias_.size(); ++n)
          if (alias_[n])
            this->bindSignal(*alias_[n], static_cast<IdxT>(n == 0 ? ReecbInternalVariables::ICMDD : ReecbInternalVariables::ICMDQ));
        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::verify() const
      {
        int  ret   = 0;
        auto check = [&](bool valid, const char* message)
        {
          if (!valid)
          {
            Log::error() << "Reecb: " << message << '\n';
            ++ret;
          }
        };
        check(S_ > ZERO<RealT> && V_ > ZERO<RealT>, "S and V must be positive");
        for (const auto value : {S_, V_, Trv_, Tp_, Vref0_, Vdip_, Vup_, dbd1_, dbd2_, kqv_, Iql1_, Iqh1_, Qmax_, Qmin_, Kqp_, Kqi_, Vmax_, Vmin_, Kvp_, Kvi_, Tiq_, Tpord_, dPmax_, dPmin_, Pmax_, Pmin_, Imax_})
          check(std::isfinite(value), "parameters must be finite");
        check(Vdip_ < Vup_, "Vdip must be less than Vup");
        check(dbd1_ <= ZERO<RealT> && ZERO<RealT> <= dbd2_, "dbd1 <= 0 <= dbd2 is required");
        check(Iql1_ <= Iqh1_, "Iql1 must not exceed Iqh1");
        check(Qmin_ <= Qmax_ && Vmin_ <= Vmax_ && Pmin_ <= Pmax_, "lower limits must not exceed upper limits");
        check(dPmin_ < ZERO<RealT> && ZERO<RealT> < dPmax_, "dPmin < 0 < dPmax is required");
        check(Imax_ > ZERO<RealT>, "Imax must be positive");
        for (const auto value : {kqv_, Kqp_, Kqi_, Kvp_, Kvi_})
          check(value >= ZERO<RealT>, "gains must be nonnegative");

        using E = ReecbExternalVariables;
        check(signals_.attachedSignals({E::VD, E::VQ, E::ID, E::IQ}).size() == 4,
              "terminal voltage and current inputs are required");
        for (size_t n = 0; n < static_cast<size_t>(E::MAXIMUM); ++n)
        {
          const auto input = static_cast<E>(n);
          for (const auto* signal : signals_.attachedSignals({input}))
            check(signal->linked(), "attached inputs must have linked sources");
        }
        return ret;
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
      {
        this->template parseInitialOutputs<Reecb>(values);
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
      {
        return this->initializeOutputs(*this, values);
      }

      template <typename scalar_type, typename index_type>
      auto Reecb<scalar_type, index_type>::operatingPoint(
          const std::map<Outputs, RealT>& outputs, const std::array<RealT, 8>& input) const -> OperatingPoint
      {
        const auto VMEAS  = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS  = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ   = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV   = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV     = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD   = static_cast<size_t>(ReecbInternalVariables::PORD);
        const auto VT     = static_cast<size_t>(ReecbInternalVariables::VT);
        const auto VSAFE  = static_cast<size_t>(ReecbInternalVariables::VSAFE);
        const auto SDIP   = static_cast<size_t>(ReecbInternalVariables::SDIP);
        const auto IQV    = static_cast<size_t>(ReecbInternalVariables::IQV);
        const auto QREF   = static_cast<size_t>(ReecbInternalVariables::QREF);
        const auto EQ     = static_cast<size_t>(ReecbInternalVariables::EQ);
        const auto VPIQ   = static_cast<size_t>(ReecbInternalVariables::VPIQ);
        const auto EPIV   = static_cast<size_t>(ReecbInternalVariables::EPIV);
        const auto RPORD  = static_cast<size_t>(ReecbInternalVariables::RPORD);
        const auto ILCAP  = static_cast<size_t>(ReecbInternalVariables::ILCAP);
        const auto IQMAX  = static_cast<size_t>(ReecbInternalVariables::IQMAX);
        const auto IPMAX  = static_cast<size_t>(ReecbInternalVariables::IPMAX);
        const auto IQBASE = static_cast<size_t>(ReecbInternalVariables::IQBASE);
        const auto IQRAW  = static_cast<size_t>(ReecbInternalVariables::IQRAW);
        const auto IQCMD  = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD  = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        this->validateOutputValues(outputs);
        for (const auto value : input)
          if (!std::isfinite(value))
            throw std::invalid_argument("Reecb: initial inputs must be finite");
        OperatingPoint point;
        point.reference = {input[4], input[5], input[6], input[7]};
        const RealT vr0 = input[0] / V_;
        const RealT vi0 = input[1] / V_;
        const RealT vt0 = std::hypot(vr0, vi0);
        if (vt0 <= ZERO<RealT>)
          throw std::invalid_argument("Reecb: initial terminal voltage must be positive");
        const RealT vmeas0      = vt0;
        const RealT vmeas_safe0 = Math::max(vmeas0, VMEAS_MINIMUM);
        const RealT idcmd0      = this->outputValue(outputs, Outputs::icmdd, input[2]);
        const RealT iqcmd_dq0   = this->outputValue(outputs, Outputs::icmdq, input[3]);
        const RealT ipcmd0      = idcmd0 * V_ / S_;
        const RealT iqcmd0      = -iqcmd_dq0 * V_ / S_;
        const RealT pmeas0      = (input[0] * input[2] + input[1] * input[3]) / S_;
        const RealT qgen0       = (input[1] * input[2] - input[0] * input[3]) / S_;
        const RealT vref0       = Vref0_given_ ? Vref0_ : vmeas0;
        if (ipcmd0 < ZERO<RealT>)
        {
          throw std::invalid_argument("Reecb: initial active-current command must be non-negative");
        }

        const RealT verr0   = Math::deadband2(vref0 - vmeas0, dbd1_, dbd2_);
        const RealT iqv0    = Math::clamp(kqv_ * verr0, Iql1_, Iqh1_);
        const RealT iqabs0  = std::abs(iqcmd0);
        RealT       iqneed0 = iqabs0;
        if (QFlag_ && iqabs0 > ZERO<RealT>)
        {
          iqneed0 += std::numbers::ln2_v<RealT> / Math::MU<RealT> + INITIALIZATION_TOLERANCE;
        }

        RealT high0 = iqabs0;
        RealT low0  = ipcmd0;
        if (Pqflag_)
        {
          high0 = ipcmd0;
          low0  = iqneed0;
        }
        // Q priority uses Imax directly for reactive current, so include the
        // smooth-clamp recovery margin carried by iqneed0.
        const auto current_limit = solveInitialLimit(
            std::max({Imax_, high0, low0, iqneed0}), high0, low0);
        if (!current_limit)
        {
          throw std::invalid_argument("Reecb: adjusted Imax cannot include the initial current commands");
        }
        const RealT imax   = current_limit->total_limit;
        const RealT ilcap0 = current_limit->off_axis_capacity;

        RealT iqmax0 = imax;
        RealT ipmax0 = ilcap0;
        if (Pqflag_)
        {
          iqmax0 = ilcap0;
          ipmax0 = imax;
        }

        RealT ipraw0 = ZERO<RealT>;
        RealT iqraw0 = ZERO<RealT>;
        if (!iclamp(ipcmd0, ZERO<RealT>, ipmax0, ipraw0)
            || !iclamp(iqcmd0, -iqmax0, iqmax0, iqraw0))
        {
          throw std::invalid_argument("Reecb: initial current commands cannot be reproduced by their limiters");
        }

        const RealT iqctl0 = iqraw0 - iqv0;
        const RealT pord0  = vmeas_safe0 * ipraw0;
        RealT       qmin   = Qmin_;
        RealT       qmax   = Qmax_;
        RealT       vmin   = Vmin_;
        RealT       vmax   = Vmax_;
        if (QFlag_ && VFlag_)
        {
          qmin = std::min(Qmin_, qgen0);
          qmax = std::max(Qmax_, qgen0);
          vmin = std::min(Vmin_, vmeas0);
          vmax = std::max(Vmax_, vmeas0);

          const RealT infinity = std::numeric_limits<RealT>::infinity();
          if (qmin == qgen0 && qmin < qmax)
          {
            qmin = std::nextafter(qmin, -infinity);
          }
          if (qmax == qgen0 && qmin < qmax)
          {
            qmax = std::nextafter(qmax, infinity);
          }
          if (vmin == vmeas0 && vmin < vmax)
          {
            vmin = std::nextafter(vmin, -infinity);
          }
          if (vmax == vmeas0 && vmin < vmax)
          {
            vmax = std::nextafter(vmax, infinity);
          }
        }
        const RealT pmin  = std::min(Pmin_, pord0);
        const RealT pmax  = std::max(Pmax_, pord0);
        const RealT pref0 = pord0 * S_;

        RealT qtarget0 = ZERO<RealT>;
        if (!QFlag_)
        {
          qtarget0 = iqctl0 * vmeas_safe0;
        }
        else if (VFlag_ && !iclamp(qgen0, qmin, qmax, qtarget0))
        {
          throw std::invalid_argument("Reecb: reactive-power limiter has no finite steady input");
        }

        RealT qref0      = ZERO<RealT>;
        RealT qext0_port = ZERO<RealT>;
        RealT pfaref0    = ZERO<RealT>;

        if (QFlag_ && !VFlag_)
        {
          qext0_port = vmeas0;
        }
        else if (PfFlag_)
        {
          if (pmeas0 == ZERO<RealT> && qtarget0 != ZERO<RealT>)
          {
            throw std::invalid_argument("Reecb: power-factor mode cannot reproduce the reactive target at zero active power");
          }
          if (pmeas0 != ZERO<RealT>)
          {
            pfaref0 = std::atan(qtarget0 / pmeas0);
          }
          qref0 = pmeas0 * std::tan(pfaref0);
          if (std::abs(qref0 - qtarget0) > std::abs(qtarget0) * INITIALIZATION_TOLERANCE)
          {
            throw std::invalid_argument("Reecb: power-factor angle cannot reproduce the reactive target");
          }
        }
        else
        {
          qext0_port = qtarget0 * S_;
          qref0      = qext0_port / S_;
        }

        const RealT eq0   = Math::clamp(qref0, qmin, qmax) - qgen0;
        RealT       xpiq0 = ZERO<RealT>;
        if (QFlag_ && VFlag_)
        {
          RealT vpiq_input0 = ZERO<RealT>;
          if (!iclamp(vmeas0, vmin, vmax, vpiq_input0))
          {
            throw std::invalid_argument("Reecb: voltage limiter has no finite steady input");
          }
          xpiq0 = vpiq_input0 - Kqp_ * eq0;
        }

        const RealT vpiq0 = Math::clamp(Kqp_ * eq0 + xpiq0, vmin, vmax);
        RealT       epiv0 = ZERO<RealT>;
        RealT       qv0   = ZERO<RealT>;
        RealT       xpiv0 = ZERO<RealT>;

        if (QFlag_)
        {
          if (VFlag_)
          {
            epiv0 = vpiq0 - vmeas0;
          }
          else
          {
            epiv0 = qext0_port - vmeas0;
          }

          if (iqmax0 <= INITIALIZATION_TOLERANCE)
          {
            xpiv0 = -Kvp_ * epiv0;
          }
          else
          {
            RealT iqctl_input0 = ZERO<RealT>;
            if (!iclamp(iqctl0, -iqmax0, iqmax0, iqctl_input0))
            {
              throw std::invalid_argument("Reecb: voltage-controller current cannot be reproduced by its limiter");
            }
            xpiv0 = iqctl_input0 - Kvp_ * epiv0;
          }
        }
        else
        {
          qv0 = qref0 / vmeas_safe0;
        }

        const RealT sdip0  = Math::inside(vt0, Vdip_, Vup_);
        RealT       qrate0 = ZERO<RealT>;
        if (QFlag_ && VFlag_)
        {
          qrate0 = sdip0 * Math::antiwindup(Kqp_ * eq0 + xpiq0, Kqi_ * eq0, vmin, vmax);
        }

        const ScalarT vstate0{Kvp_ * epiv0 + xpiv0};
        const ScalarT vderiv0{Kvi_ * epiv0};
        RealT         vrate0 = ZERO<RealT>;
        if (QFlag_)
        {
          vrate0 = sdip0 * static_cast<RealT>(awband(vstate0, vderiv0, ScalarT{iqmax0}));
        }

        const RealT iqbase0     = Math::clamp(Kvp_ * epiv0 + xpiv0, -iqmax0, iqmax0);
        RealT       iqraw_check = qv0 + iqv0;
        if (QFlag_)
        {
          iqraw_check = iqbase0 + iqv0;
        }

        const RealT iqcmd_check = Math::clamp(iqraw_check, -iqmax0, iqmax0);
        const RealT ipcmd_check = Math::clamp(pord0 / vmeas_safe0, ZERO<RealT>, ipmax0);

        if (!std::isfinite(imax) || !std::isfinite(ilcap0)
            || !std::isfinite(iqmax0) || !std::isfinite(ipmax0)
            || !std::isfinite(ipraw0) || !std::isfinite(iqraw0) || !std::isfinite(pord0)
            || !std::isfinite(pref0) || !std::isfinite(qtarget0) || !std::isfinite(qref0)
            || !std::isfinite(qext0_port) || !std::isfinite(pfaref0) || !std::isfinite(eq0)
            || !std::isfinite(xpiq0) || !std::isfinite(epiv0) || !std::isfinite(xpiv0)
            || !std::isfinite(qv0) || !std::isfinite(qrate0) || !std::isfinite(vrate0)
            || !std::isfinite(iqbase0) || !std::isfinite(iqraw_check)
            || !std::isfinite(iqcmd_check) || !std::isfinite(ipcmd_check))
        {
          throw std::invalid_argument("Reecb: initialization produced a nonfinite value");
        }
        if (std::abs(qrate0) > INITIALIZATION_TOLERANCE || std::abs(vrate0) > INITIALIZATION_TOLERANCE)
        {
          throw std::invalid_argument("Reecb: controller state rate is nonzero at initialization");
        }
        if (std::abs(iqcmd_check - iqcmd0) > INITIALIZATION_TOLERANCE
            || std::abs(ipcmd_check - ipcmd0) > INITIALIZATION_TOLERANCE)
        {
          throw std::invalid_argument("Reecb: current-command limiter reconstruction is inexact");
        }

        point.state[VMEAS]  = vmeas0;
        point.state[PMEAS]  = pmeas0;
        point.state[XPIQ]   = xpiq0;
        point.state[XPIV]   = xpiv0;
        point.state[QV]     = qv0;
        point.state[PORD]   = pord0;
        point.state[VT]     = vt0;
        point.state[VSAFE]  = vmeas_safe0;
        point.state[SDIP]   = sdip0;
        point.state[IQV]    = iqv0;
        point.state[QREF]   = qref0;
        point.state[EQ]     = eq0;
        point.state[VPIQ]   = vpiq0;
        point.state[EPIV]   = epiv0;
        point.state[RPORD]  = ZERO<RealT>;
        point.state[ILCAP]  = ilcap0;
        point.state[IQMAX]  = iqmax0;
        point.state[IPMAX]  = ipmax0;
        point.state[IQBASE] = iqbase0;
        point.state[IQRAW]  = iqraw_check;
        point.state[IQCMD]  = iqcmd0;
        point.state[IPCMD]  = ipcmd0;
        point.state[22]     = idcmd0;
        point.state[23]     = iqcmd_dq0;
        point.reference[0]  = pref0;
        if (QFlag_ && !VFlag_)
          point.reference[2] = qext0_port * V_;
        else if (PfFlag_)
          point.reference[3] = pfaref0;
        else
          point.reference[1] = qext0_port;

        point.qmin = qmin;
        point.qmax = qmax;
        point.vmin = vmin;
        point.vmax = vmax;
        point.pmin = pmin;
        point.pmax = pmax;
        point.imax = imax;
        point.vref = vref0;
        return point;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        if (const int errors = verify(); errors != 0)
          return errors;
        using E = ReecbExternalVariables;
        std::array<RealT, 8> input{};
        for (size_t n = 0; n < input.size(); ++n)
          for (const auto* signal : signals_.attachedSignals({static_cast<E>(n)}))
            input[n] = static_cast<RealT>(signal->read());
        const auto point = operatingPoint(outputs, input);
        for (size_t n = 0; n < point.reference.size(); ++n)
          if (referenceActive(n) && !signals_.attachedSignals({static_cast<E>(4 + n)}).empty()
              && std::abs(input[4 + n] - point.reference[n]) > RealT{1e-10} * (ONE<RealT> + std::abs(point.reference[n])))
            throw std::invalid_argument("Reecb: initial reference is inconsistent with the requested output");
        reference_ = point.reference;
        Vmin_      = point.vmin;
        Vmax_      = point.vmax;
        Imax_      = point.imax;
        Vref0_     = point.vref;
        Qmin_      = point.qmin;
        Qmax_      = point.qmax;
        Pmin_      = point.pmin;
        Pmax_      = point.pmax;
        for (size_t n = 0; n < point.state.size(); ++n)
          y_.getData()[n] = static_cast<ScalarT>(point.state[n]);
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        y_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      typename Component<scalar_type, index_type>::InitializationPortsT Reecb<scalar_type, index_type>::initializationPorts()
      {
        using E = ReecbExternalVariables;
        typename Component<ScalarT, IdxT>::InitializationPortsT ports;
        ports.inputs = signals_.attachedSignals({E::VD, E::VQ, E::ID, E::IQ, E::PREF, E::QREF, E::VREF, E::PFAREF});
        for (size_t n = 0; n < reference_.size(); ++n)
          if (referenceActive(n))
            for (auto* signal : signals_.attachedSignals({static_cast<E>(4 + n)}))
              ports.targets.push_back(signal);
        for (size_t n = 0; n < output_.size(); ++n)
        {
          const auto name = std::string(magic_enum::enum_name(static_cast<Outputs>(n)));
          ports.outputs.emplace(name, &output_[n]);
          if (alias_[n])
            ports.outputs.emplace(name, alias_[n]);
        }
        return ports;
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::prepareInitialization(typename Component<ScalarT, IdxT>::InitialStateT& initial)
      {
        using E = ReecbExternalVariables;
        std::array<RealT, 8> input{};
        for (size_t n = 0; n < 4; ++n)
          for (const auto* signal : signals_.attachedSignals({static_cast<E>(n)}))
            input[n] = initial.value(*signal);
        const auto outputs = this->template parseInitialOutputs<Reecb>(initial.outputs(*this));
        const auto point   = operatingPoint(outputs, input);
        for (size_t n = 0; n < point.reference.size(); ++n)
          if (referenceActive(n))
            for (const auto* signal : signals_.attachedSignals({static_cast<E>(4 + n)}))
              initial.require(*signal, point.reference[n], *this);
        initial.provide(output_[0], point.state[22]);
        initial.provide(output_[1], point.state[23]);
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        abs_tol_.getData()[22] = static_cast<ScalarT>(tolerance * S_ / V_);
        abs_tol_.getData()[23] = static_cast<ScalarT>(tolerance * S_ / V_);
        abs_tol_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline int
      Reecb<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* input,
          const ScalarT*,
          ScalarT* f)
      {
        const auto VMEAS  = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS  = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ   = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV   = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV     = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD   = static_cast<size_t>(ReecbInternalVariables::PORD);
        const auto VT     = static_cast<size_t>(ReecbInternalVariables::VT);
        const auto VSAFE  = static_cast<size_t>(ReecbInternalVariables::VSAFE);
        const auto SDIP   = static_cast<size_t>(ReecbInternalVariables::SDIP);
        const auto IQV    = static_cast<size_t>(ReecbInternalVariables::IQV);
        const auto QREF   = static_cast<size_t>(ReecbInternalVariables::QREF);
        const auto EQ     = static_cast<size_t>(ReecbInternalVariables::EQ);
        const auto VPIQ   = static_cast<size_t>(ReecbInternalVariables::VPIQ);
        const auto EPIV   = static_cast<size_t>(ReecbInternalVariables::EPIV);
        const auto RPORD  = static_cast<size_t>(ReecbInternalVariables::RPORD);
        const auto ILCAP  = static_cast<size_t>(ReecbInternalVariables::ILCAP);
        const auto IQMAX  = static_cast<size_t>(ReecbInternalVariables::IQMAX);
        const auto IPMAX  = static_cast<size_t>(ReecbInternalVariables::IPMAX);
        const auto IQBASE = static_cast<size_t>(ReecbInternalVariables::IQBASE);
        const auto IQRAW  = static_cast<size_t>(ReecbInternalVariables::IQRAW);
        const auto IQCMD  = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD  = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        const ScalarT vmeas  = y[VMEAS];
        const ScalarT pmeas  = y[PMEAS];
        const ScalarT xpiq   = y[XPIQ];
        const ScalarT xpiv   = y[XPIV];
        const ScalarT qv     = y[QV];
        const ScalarT pord   = y[PORD];
        const ScalarT vt     = y[VT];
        const ScalarT vsafe  = y[VSAFE];
        const ScalarT sdip   = y[SDIP];
        const ScalarT iqv    = y[IQV];
        const ScalarT qref   = y[QREF];
        const ScalarT eq     = y[EQ];
        const ScalarT vpiq   = y[VPIQ];
        const ScalarT epiv   = y[EPIV];
        const ScalarT rpord  = y[RPORD];
        const ScalarT ilcap  = y[ILCAP];
        const ScalarT iqmax  = y[IQMAX];
        const ScalarT ipmax  = y[IPMAX];
        const ScalarT iqbase = y[IQBASE];
        const ScalarT iqraw  = y[IQRAW];
        const ScalarT iqcmd  = y[IQCMD];
        const ScalarT ipcmd  = y[IPCMD];

        const ScalarT vmeas_dot = yp[VMEAS];
        const ScalarT pmeas_dot = yp[PMEAS];
        const ScalarT xpiq_dot  = yp[XPIQ];
        const ScalarT xpiv_dot  = yp[XPIV];
        const ScalarT qv_dot    = yp[QV];
        const ScalarT pord_dot  = yp[PORD];

        const ScalarT vr     = input[0] / V_;
        const ScalarT vi     = input[1] / V_;
        const ScalarT pe     = (input[0] * input[2] + input[1] * input[3]) / S_;
        const ScalarT qgen   = (input[1] * input[2] - input[0] * input[3]) / S_;
        const ScalarT pref   = input[4] / S_;
        const ScalarT qext   = input[5] / S_;
        const ScalarT vref   = input[6] / V_;
        const ScalarT pfaref = input[7];

        const ScalarT verr        = Math::deadband2(Vref0_ - vmeas, dbd1_, dbd2_);
        const ScalarT q_pi_state  = Kqp_ * eq + xpiq;
        const ScalarT v_pi_state  = Kvp_ * epiv + xpiv;
        const ScalarT fpord       = (pref - pord) / Tpord_;
        // Select before the factored square to avoid 0 * inf on the inactive path.
        const ScalarT high        = pq_on_ * ipcmd + pq_off_ * iqcmd;
        const ScalarT q_pi_rate   = q_pi_on_ * sdip * Math::antiwindup(q_pi_state, Kqi_ * eq, Vmin_, Vmax_);
        const ScalarT v_pi_rate   = q_on_ * sdip * awband(v_pi_state, Kvi_ * epiv, iqmax);
        const ScalarT qv_rate     = q_off_ * sdip * (qref / vsafe - qv) / Tiq_;
        const ScalarT pord_rate   = sdip * Math::antiwindup(pord, rpord, Pmin_, Pmax_);
        const ScalarT iqv_target  = Math::clamp(kqv_ * verr, Iql1_, Iqh1_);
        const ScalarT qref_target = q_ref_on_ * (pf_on_ * pmeas * std::tan(pfaref) + pf_off_ * qext);

        f[VMEAS]  = -vmeas_dot + (vt - vmeas) / Trv_;
        f[PMEAS]  = -pmeas_dot + (pe - pmeas) / Tp_;
        f[XPIQ]   = -xpiq_dot + q_pi_rate;
        f[XPIV]   = -xpiv_dot + v_pi_rate;
        f[QV]     = -qv_dot + qv_rate;
        f[PORD]   = -pord_dot + pord_rate;
        f[VT]     = -vt * vt + vr * vr + vi * vi;
        f[VSAFE]  = -vsafe + Math::max(vmeas, VMEAS_MINIMUM);
        f[SDIP]   = -sdip + Math::inside(vt, Vdip_, Vup_);
        f[IQV]    = -iqv + iqv_target;
        f[QREF]   = -qref + qref_target;
        f[EQ]     = -eq + Math::clamp(qref, Qmin_, Qmax_) - qgen;
        f[VPIQ]   = -vpiq + Math::clamp(q_pi_state, Vmin_, Vmax_);
        f[EPIV]   = -epiv + q_pi_on_ * vpiq + v_ref_on_ * vref - q_on_ * vmeas;
        f[RPORD]  = -rpord + aslew(fpord, dPmin_, dPmax_);
        f[ILCAP]  = -ilcap + sqrtramp(circleSquare(Imax_, high));
        f[IQMAX]  = -iqmax + pq_on_ * ilcap + pq_off_ * Imax_;
        f[IPMAX]  = -ipmax + pq_on_ * Imax_ + pq_off_ * ilcap;
        f[IQBASE] = -iqbase + Math::clamp(v_pi_state, -iqmax, iqmax);
        f[IQRAW]  = -iqraw + q_on_ * iqbase + q_off_ * qv + iqv;
        f[IQCMD]  = -iqcmd + Math::clamp(iqraw, -iqmax, iqmax);
        f[IPCMD]  = -ipcmd + Math::clamp(pord / vsafe, ZERO<RealT>, ipmax);

        f[22] = -y[22] + (S_ / V_) * ipcmd;
        f[23] = -y[23] - (S_ / V_) * iqcmd;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateInternalResidual()
      {
        for (size_t n = 0; n < reference_.size(); ++n)
          y_ext_[4 + n] = static_cast<ScalarT>(reference_[n]);
        this->gatherExternalVariables();
        evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initializeMonitor()
      {
        using Mon = typename ModelDataT::MonitorableVariables;
        monitor_->set(Mon::icmdd, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::ICMDD)]; });
        monitor_->set(Mon::icmdq, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::ICMDQ)]; });
        monitor_->set(Mon::ipcmd, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::IPCMD)]; });
        monitor_->set(Mon::iqcmd, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::IQCMD)]; });
        monitor_->set(Mon::iqv, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::IQV)]; });
        monitor_->set(Mon::vmeas, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::VMEAS)]; });
        monitor_->set(Mon::pmeas, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::PMEAS)]; });
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Reecb<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      bool Reecb<scalar_type, index_type>::referenceActive(size_t n) const
      {
        return n == 0 || (n == 1 && !(QFlag_ && !VFlag_) && !PfFlag_)
               || (n == 2 && QFlag_ && !VFlag_) || (n == 3 && !(QFlag_ && !VFlag_) && PfFlag_);
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      Reecb<scalar_type, index_type>::aslew(ScalarT rate, RealT lower, RealT upper)
      {
        assert(lower < ZERO<RealT> && ZERO<RealT> < upper);
        return rate
               / (ONE<RealT> + Math::ramp(rate / upper - ONE<RealT>) + Math::ramp(rate / lower - ONE<RealT>));
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      Reecb<scalar_type, index_type>::awband(ScalarT state, ScalarT rate, ScalarT band)
      {
        const ScalarT above_min = Math::above(state, -band);
        const ScalarT below_max = Math::below(state, band);
        return (above_min * below_max + (ONE<RealT> - below_max) * Math::sigmoid(-rate)
                + (ONE<RealT> - above_min) * Math::sigmoid(rate))
               * rate;
      }

      template <typename scalar_type, typename index_type>
      template <typename ValueT>
      [[gnu::always_inline]] inline ValueT
      Reecb<scalar_type, index_type>::sqrtramp(ValueT x)
      {
        const RealT root_width = ONE<RealT> / Math::MU<RealT>;
        // Keep closed-circle leakage below the initialization tolerance.
        const RealT knee       = INITIALIZATION_TOLERANCE / Math::MU<RealT>;

        const ValueT absolute    = std::abs(x);
        const ValueT normalizer  = absolute + knee;
        const ValueT scaled_x    = absolute / normalizer;
        const ValueT scaled_knee = knee / normalizer;
        const ValueT magnitude   = normalizer
                                 * std::sqrt(scaled_x * scaled_x
                                             + scaled_knee * scaled_knee);
        const ValueT half_sum = HALF<RealT> * magnitude
                                + HALF<RealT> * absolute;
        const ValueT conjugate = QUARTER<RealT> * knee * (knee / half_sum);
        const ValueT hinged    = HALF<RealT> * x
                              + HALF<RealT> * absolute
                              + conjugate;
        // Divide through by MU to avoid overflow for large finite circles.
        return hinged / std::sqrt(hinged + root_width * root_width);
      }

      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::isqrtramp(RealT y)
      {
        if (y <= sqrtramp(ZERO<RealT>))
        {
          return ZERO<RealT>;
        }

        const RealT mu       = Math::MU<RealT>;
        const RealT knee     = INITIALIZATION_TOLERANCE / mu;
        const RealT scaled_y = HALF<RealT> * mu * y;
        const RealT hinged   = y / mu
                             * (scaled_y + std::hypot(scaled_y, ONE<RealT>));
        const RealT square = hinged - QUARTER<RealT> * knee * (knee / hinged);
        return std::max(ZERO<RealT>, square);
      }

      template <typename scalar_type, typename index_type>
      template <typename ValueT>
      [[gnu::always_inline]] inline ValueT
      Reecb<scalar_type, index_type>::circleSquare(RealT limit, ValueT high)
      {
        return ((HALF<RealT> * limit - HALF<RealT> * high)
                * (HALF<RealT> * limit + HALF<RealT> * high))
               * FOUR<RealT>;
      }

      template <typename scalar_type, typename index_type>
      auto Reecb<scalar_type, index_type>::solveInitialLimit(
          RealT lower, RealT high, RealT low) -> std::optional<InitialCurrentLimit>
      {
        RealT imax = lower;
        RealT cap  = sqrtramp(circleSquare(imax, high));
        if (!std::isfinite(cap))
        {
          return std::nullopt;
        }
        if (cap >= low)
        {
          return InitialCurrentLimit{imax, cap};
        }

        const RealT square = isqrtramp(low);
        if (!std::isfinite(square))
        {
          return std::nullopt;
        }

        imax = std::max(lower, std::hypot(high, std::sqrt(square)));
        cap  = sqrtramp(circleSquare(imax, high));

        for (int step = 0;
             step < std::numeric_limits<RealT>::digits && cap < low;
             ++step)
        {
          const RealT next = std::nextafter(imax, std::numeric_limits<RealT>::max());
          if (!std::isfinite(next) || next <= imax)
          {
            return std::nullopt;
          }
          imax = next;
          cap  = sqrtramp(circleSquare(imax, high));
          if (!std::isfinite(cap))
          {
            return std::nullopt;
          }
        }

        if (!std::isfinite(imax) || !std::isfinite(cap) || cap < low)
        {
          return std::nullopt;
        }
        return InitialCurrentLimit{imax, cap};
      }

      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::logOneMinusExp(RealT x)
      {
        static constexpr auto log_two = std::numbers::ln2_v<RealT>;

        if (x < log_two)
        {
          return log_two - HALF<RealT> * x + std::log(std::sinh(HALF<RealT> * x));
        }
        return std::log1p(-std::exp(-x));
      }

      template <typename scalar_type, typename index_type>
      bool Reecb<scalar_type, index_type>::iclamp(RealT output, RealT lower, RealT upper, RealT& input) const
      {
        if (!std::isfinite(output) || !std::isfinite(lower) || !std::isfinite(upper) || lower > upper
            || output < lower - INITIALIZATION_TOLERANCE || output > upper + INITIALIZATION_TOLERANCE)
        {
          return false;
        }

        output = std::clamp(output, lower, upper);
        if (upper == lower)
        {
          input = lower;
          return true;
        }

        const RealT mu     = Math::MU<RealT>;
        const RealT offset = -std::log(std::expm1(mu * HALF<RealT> * INITIALIZATION_TOLERANCE)) / mu;
        if (output == lower)
        {
          input = lower - offset;
          return true;
        }
        if (output == upper)
        {
          input = upper + offset;
          return true;
        }

        const RealT a = mu * (output - lower);
        const RealT b = mu * (upper - output);
        input         = lower + (a + logOneMinusExp(a) - logOneMinusExp(b)) / mu;
        return std::isfinite(input);
      }

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

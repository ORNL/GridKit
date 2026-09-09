#pragma once

#include <cmath>
#include <mutex>
#include <numbers>

#include <GridKit/Model/EMT/Component/Controller/REPCA/Repca.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca()
        : Repca(ModelDataT{})
      {
      }

      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        size_ = static_cast<IdxT>(RepcaInternalVariables::MAXIMUM);
        signals_.template assignSignal<RepcaInternalVariables::QEXT>(&output_[0]);
        signals_.template assignSignal<RepcaInternalVariables::PEXT>(&output_[1]);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::~Repca() = default;

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;
        S_              = parameter<RealT>(data, Parameter::S, S_);
        V_              = parameter<RealT>(data, Parameter::V, V_);
        VcompFlag_      = parameter<bool>(data, Parameter::VcompFlag, VcompFlag_);
        RefFlag_        = parameter<bool>(data, Parameter::RefFlag, RefFlag_);
        Freqflag_       = parameter<bool>(data, Parameter::Freqflag, Freqflag_);
        Tfltr_          = parameter<RealT>(data, Parameter::Tfltr, Tfltr_);
        Vfrz_           = parameter<RealT>(data, Parameter::Vfrz, Vfrz_);
        Rc_             = parameter<RealT>(data, Parameter::Rc, Rc_);
        Xc_             = parameter<RealT>(data, Parameter::Xc, Xc_);
        Kc_             = parameter<RealT>(data, Parameter::Kc, Kc_);
        dbdlow_         = parameter<RealT>(data, Parameter::dbdlow, dbdlow_);
        dbdupper_       = parameter<RealT>(data, Parameter::dbdupper, dbdupper_);
        emax_           = parameter<RealT>(data, Parameter::emax, emax_);
        emin_           = parameter<RealT>(data, Parameter::emin, emin_);
        Kp_             = parameter<RealT>(data, Parameter::Kp, Kp_);
        Ki_             = parameter<RealT>(data, Parameter::Ki, Ki_);
        Qmax_           = parameter<RealT>(data, Parameter::Qmax, Qmax_);
        Qmin_           = parameter<RealT>(data, Parameter::Qmin, Qmin_);
        Tft_            = parameter<RealT>(data, Parameter::Tft, Tft_);
        Tfv_            = parameter<RealT>(data, Parameter::Tfv, Tfv_);
        Tp_             = parameter<RealT>(data, Parameter::Tp, Tp_);
        fdbd1_          = parameter<RealT>(data, Parameter::fdbd1, fdbd1_);
        fdbd2_          = parameter<RealT>(data, Parameter::fdbd2, fdbd2_);
        Ddn_            = parameter<RealT>(data, Parameter::Ddn, Ddn_);
        Dup_            = parameter<RealT>(data, Parameter::Dup, Dup_);
        femax_          = parameter<RealT>(data, Parameter::femax, femax_);
        femin_          = parameter<RealT>(data, Parameter::femin, femin_);
        Kpg_            = parameter<RealT>(data, Parameter::Kpg, Kpg_);
        Kig_            = parameter<RealT>(data, Parameter::Kig, Kig_);
        Pmax_           = parameter<RealT>(data, Parameter::Pmax, Pmax_);
        Pmin_           = parameter<RealT>(data, Parameter::Pmin, Pmin_);
        Tlag_           = parameter<RealT>(data, Parameter::Tlag, Tlag_);
        for (const auto value : {Tfltr_, Tft_, Tfv_, Tp_, Tlag_})
          if (!std::isfinite(value) || value < ZERO<RealT>)
            throw std::invalid_argument("Repca: time constants must be finite and nonnegative");
        if (std::min({Tfltr_, Tfv_, Tp_, Tlag_}) < TIME_CONSTANT_MINIMUM)
        {
          static std::once_flag warning;
          std::call_once(warning, []
                         { Log::warning() << "Repca: Tfltr, Tfv, Tp, and Tlag below 0.001 s are raised to that floor\n"; });
        }
        Tfltr_     = std::max(Tfltr_, TIME_CONSTANT_MINIMUM);
        Tfv_       = std::max(Tfv_, TIME_CONSTANT_MINIMUM);
        Tp_        = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        Tlag_      = std::max(Tlag_, TIME_CONSTANT_MINIMUM);
        vcomp_on_  = static_cast<RealT>(VcompFlag_);
        vcomp_off_ = ONE<RealT> - vcomp_on_;
        ref_on_    = static_cast<RealT>(RefFlag_);
        ref_off_   = ONE<RealT> - ref_on_;
        freq_on_   = static_cast<RealT>(Freqflag_);
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::attachInput(InputSignals inputs)
      {
        if (allocated_)
          throw std::logic_error("Repca inputs cannot change after allocation");
        for (size_t p = 0; p < inputs.size(); ++p)
          if (inputs[p])
            signals_.attachSignal(static_cast<RepcaExternalVariables>(p), inputs[p]);
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        if (allocated_)
          throw std::logic_error("Repca outputs cannot change after allocation");
        const auto n = static_cast<size_t>(output);
        if (n >= output_.size() || !signal || alias_[n])
          throw std::invalid_argument("Repca: invalid or duplicate output assignment");
        signal->claimProducer();
        alias_[n] = signal;
      }

      template <typename scalar_type, typename index_type>
      typename Repca<scalar_type, index_type>::SignalT& Repca<scalar_type, index_type>::outputSignal(Outputs output)
      {
        return output_.at(static_cast<size_t>(output));
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
          this->allocateVectors(size_);
        tag_.resize(static_cast<size_t>(size_));
        variable_indices_.resize(static_cast<size_t>(size_));
        residual_indices_.resize(static_cast<size_t>(size_));
        this->assignGlobalIndices(0);
        this->allocateExternalVectors(static_cast<IdxT>(RepcaExternalVariables::MAXIMUM), 0);
        signals_.registerExternalVariableSignals(*this);
        signals_.bindInternalVariableSignals(*this);
        for (size_t n = 0; n < alias_.size(); ++n)
          if (alias_[n])
            this->bindSignal(*alias_[n], static_cast<IdxT>(n == 0 ? RepcaInternalVariables::QEXT : RepcaInternalVariables::PEXT));
        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::verify() const
      {
        int  ret   = 0;
        auto check = [&](bool valid, const char* message)
        {
          if (!valid)
          {
            Log::error() << "Repca: " << message << '\n';
            ++ret;
          }
        };
        check(S_ > ZERO<RealT> && V_ > ZERO<RealT>, "S and V must be positive");
        for (const auto value : {S_, V_, Tfltr_, Vfrz_, Rc_, Xc_, Kc_, dbdlow_, dbdupper_, emax_, emin_, Kp_, Ki_, Qmax_, Qmin_, Tft_, Tfv_, Tp_, fdbd1_, fdbd2_, Ddn_, Dup_, femax_, femin_, Kpg_, Kig_, Pmax_, Pmin_, Tlag_, vcomp_on_, vcomp_off_, ref_on_, ref_off_, freq_on_})
          check(std::isfinite(value), "parameters must be finite");
        check(dbdlow_ <= ZERO<RealT> && ZERO<RealT> <= dbdupper_,
              "dbdlow <= 0 <= dbdupper is required");
        check(emin_ <= ZERO<RealT> && ZERO<RealT> <= emax_,
              "emin <= 0 <= emax is required");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(fdbd1_ <= ZERO<RealT> && ZERO<RealT> <= fdbd2_,
              "fdbd1 <= 0 <= fdbd2 is required");
        check(Ddn_ >= ZERO<RealT>, "Ddn must be non-negative");
        check(Dup_ >= ZERO<RealT>, "Dup must be non-negative");
        check(femin_ <= ZERO<RealT> && ZERO<RealT> <= femax_,
              "femin <= 0 <= femax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");

        using E = RepcaExternalVariables;
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
      void Repca<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
      {
        this->template parseInitialOutputs<Repca>(values);
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
      {
        return this->initializeOutputs(*this, values);
      }

      template <typename scalar_type, typename index_type>
      auto Repca<scalar_type, index_type>::operatingPoint(
          const std::map<Outputs, RealT>& outputs, const std::array<RealT, 9>& input) const -> OperatingPoint
      {
        this->validateOutputValues(outputs);
        for (const RealT value : input)
          if (!std::isfinite(value))
            throw std::invalid_argument("Repca: initial inputs must be finite");
        const RealT vd     = input[0] / V_;
        const RealT vq     = input[1] / V_;
        const RealT id     = input[2] * V_ / S_;
        const RealT iq     = input[3] * V_ / S_;
        const RealT p      = vd * id + vq * iq;
        const RealT q      = vq * id - vd * iq;
        const RealT v      = std::hypot(vd, vq);
        const RealT vldc   = std::hypot(vd - Rc_ * id + Xc_ * iq, vq - Rc_ * iq - Xc_ * id);
        const RealT vdroop = v + Kc_ * q;
        const RealT vctrl  = vcomp_on_ * vldc + vcomp_off_ * vdroop;
        const RealT qext   = this->outputValue(outputs, Outputs::qext, q * S_);
        const RealT pext   = this->outputValue(outputs, Outputs::pext, Freqflag_ ? p * S_ : ZERO<RealT>);
        if (!Freqflag_ && pext != ZERO<RealT>)
          throw std::invalid_argument("Repca: pext must be zero when Freqflag is false");
        const RealT qpi  = qext / S_;
        const RealT ppi  = Freqflag_ ? pext / S_ : p;
        const RealT qmin = std::min(Qmin_, qpi);
        const RealT qmax = std::max(Qmax_, qpi);
        const RealT pmin = std::min(Pmin_, ppi);
        const RealT pmax = std::max(Pmax_, ppi);
        RealT       erqdb, erq, ep, efref, xqpi, xppi;
        if (!invertClamp(ZERO<RealT>, emin_, emax_, erqdb)
            || !invertDeadband(erqdb, dbdlow_, dbdupper_, erq)
            || !invertClamp(ZERO<RealT>, femin_, femax_, ep)
            || !invertDeadband(ZERO<RealT>, fdbd1_, fdbd2_, efref)
            || !invertClamp(qpi, qmin, qmax, xqpi)
            || !invertClamp(ppi, pmin, pmax, xppi))
          throw std::invalid_argument("Repca: limiters have no finite steady input");
        OperatingPoint point{
            {vctrl, q, xqpi, qpi, p, xppi, ppi, v, vldc, vdroop, vctrl, Math::above(v, Vfrz_), erq, erqdb, 0, qpi, qext, 0, ep, 0, ppi, pext},
            {(vctrl + (RefFlag_ ? erq : ZERO<RealT>) ) * V_, (p + ep) * S_, (q + (RefFlag_ ? ZERO<RealT> : erq)) * S_, input[4] + efref},
            qmin,
            qmax,
            pmin,
            pmax};
        for (const RealT value : point.state)
          if (!std::isfinite(value))
            throw std::invalid_argument("Repca: nonfinite derived initial state");
        for (const RealT value : point.reference)
          if (!std::isfinite(value))
            throw std::invalid_argument("Repca: nonfinite derived initial reference");
        return point;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        if (const int errors = verify(); errors != 0)
          return errors;
        using E = RepcaExternalVariables;
        std::array<RealT, 9> input{0, 0, 0, 0, 1, 0, 0, 0, 1};
        for (size_t n = 0; n < input.size(); ++n)
          for (const auto* signal : signals_.attachedSignals({static_cast<E>(n)}))
            input[n] = static_cast<RealT>(signal->read());
        const auto point = operatingPoint(outputs, input);
        for (size_t n = 0; n < point.reference.size(); ++n)
          if (!signals_.attachedSignals({static_cast<E>(5 + n)}).empty()
              && std::abs(input[5 + n] - point.reference[n]) > RealT{1e-10} * (ONE<RealT> + std::abs(point.reference[n])))
            throw std::invalid_argument("Repca: initial reference is inconsistent with the requested output");
        reference_ = point.reference;
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
      typename Component<scalar_type, index_type>::InitializationPortsT Repca<scalar_type, index_type>::initializationPorts()
      {
        using E = RepcaExternalVariables;
        typename Component<ScalarT, IdxT>::InitializationPortsT ports;
        ports.inputs  = signals_.attachedSignals({E::VD, E::VQ, E::ID, E::IQ, E::FREQ, E::VREF, E::PREF, E::QREF, E::FREQREF});
        ports.targets = signals_.attachedSignals({E::VREF, E::PREF, E::QREF, E::FREQREF});
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
      void Repca<scalar_type, index_type>::prepareInitialization(typename Component<ScalarT, IdxT>::InitialStateT& initial)
      {
        using E = RepcaExternalVariables;
        std::array<RealT, 9> input{0, 0, 0, 0, 1, 0, 0, 0, 1};
        for (size_t n = 0; n < 5; ++n)
          for (const auto* signal : signals_.attachedSignals({static_cast<E>(n)}))
            input[n] = initial.value(*signal);
        const auto outputs = this->template parseInitialOutputs<Repca>(initial.outputs(*this));
        const auto point   = operatingPoint(outputs, input);
        for (size_t n = 0; n < point.reference.size(); ++n)
          for (const auto* signal : signals_.attachedSignals({static_cast<E>(5 + n)}))
            initial.require(*signal, point.reference[n], *this);
        initial.provide(output_[0], point.state[16]);
        initial.provide(output_[1], point.state[21]);
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        abs_tol_.getData()[16] = static_cast<ScalarT>(tolerance * S_);
        abs_tol_.getData()[21] = static_cast<ScalarT>(tolerance * S_);
        abs_tol_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y, const ScalarT* yp, const ScalarT* input, const ScalarT*, ScalarT* f)
      {
        const auto VMEAS      = static_cast<size_t>(RepcaInternalVariables::VMEAS);
        const auto QMEAS      = static_cast<size_t>(RepcaInternalVariables::QMEAS);
        const auto XQPI       = static_cast<size_t>(RepcaInternalVariables::XQPI);
        const auto XQLAG      = static_cast<size_t>(RepcaInternalVariables::XQLAG);
        const auto PMEAS      = static_cast<size_t>(RepcaInternalVariables::PMEAS);
        const auto XPPI       = static_cast<size_t>(RepcaInternalVariables::XPPI);
        const auto PREF_STATE = static_cast<size_t>(RepcaInternalVariables::PREF);
        const auto V          = static_cast<size_t>(RepcaInternalVariables::V);
        const auto VLDC       = static_cast<size_t>(RepcaInternalVariables::VLDC);
        const auto VDROOP     = static_cast<size_t>(RepcaInternalVariables::VDROOP);
        const auto VCTRL      = static_cast<size_t>(RepcaInternalVariables::VCTRL);
        const auto SFRZ       = static_cast<size_t>(RepcaInternalVariables::SFRZ);
        const auto ERQ        = static_cast<size_t>(RepcaInternalVariables::ERQ);
        const auto ERQDB      = static_cast<size_t>(RepcaInternalVariables::ERQDB);
        const auto ERQLIM     = static_cast<size_t>(RepcaInternalVariables::ERQLIM);
        const auto QPI        = static_cast<size_t>(RepcaInternalVariables::QPI);
        const auto QEXT       = static_cast<size_t>(RepcaInternalVariables::QEXT);
        const auto EF         = static_cast<size_t>(RepcaInternalVariables::EF);
        const auto EP         = static_cast<size_t>(RepcaInternalVariables::EP);
        const auto EPLIM      = static_cast<size_t>(RepcaInternalVariables::EPLIM);
        const auto PPI        = static_cast<size_t>(RepcaInternalVariables::PPI);
        const auto PEXT       = static_cast<size_t>(RepcaInternalVariables::PEXT);

        const ScalarT vmeas  = y[VMEAS];
        const ScalarT qmeas  = y[QMEAS];
        const ScalarT xqpi   = y[XQPI];
        const ScalarT xqlag  = y[XQLAG];
        const ScalarT pmeas  = y[PMEAS];
        const ScalarT xppi   = y[XPPI];
        const ScalarT pref   = y[PREF_STATE];
        const ScalarT v      = y[V];
        const ScalarT vldc   = y[VLDC];
        const ScalarT vdroop = y[VDROOP];
        const ScalarT vctrl  = y[VCTRL];
        const ScalarT sfrz   = y[SFRZ];
        const ScalarT erq    = y[ERQ];
        const ScalarT erqdb  = y[ERQDB];
        const ScalarT erqlim = y[ERQLIM];
        const ScalarT qpi    = y[QPI];
        const ScalarT qext   = y[QEXT] / S_;
        const ScalarT ef     = y[EF];
        const ScalarT ep     = y[EP];
        const ScalarT eplim  = y[EPLIM];
        const ScalarT ppi    = y[PPI];
        const ScalarT pext   = y[PEXT] / S_;

        const ScalarT vmeas_dot = yp[VMEAS];
        const ScalarT qmeas_dot = yp[QMEAS];
        const ScalarT xqpi_dot  = yp[XQPI];
        const ScalarT xqlag_dot = yp[XQLAG];
        const ScalarT pmeas_dot = yp[PMEAS];
        const ScalarT xppi_dot  = yp[XPPI];
        const ScalarT pref_dot  = yp[PREF_STATE];

        const ScalarT vr      = input[0] / V_;
        const ScalarT vi      = input[1] / V_;
        const ScalarT ir      = input[2] * V_ / S_;
        const ScalarT ii      = input[3] * V_ / S_;
        const ScalarT p       = vr * ir + vi * ii;
        const ScalarT q       = vi * ir - vr * ii;
        const ScalarT freq    = input[4];
        const ScalarT vref    = input[5] / V_;
        const ScalarT pref_in = input[6] / S_;
        const ScalarT qref    = input[7] / S_;
        const ScalarT freqref = input[8];

        const ScalarT vldc_r = vr - Rc_ * ir + Xc_ * ii;
        const ScalarT vldc_i = vi - Rc_ * ii - Xc_ * ir;
        const ScalarT pfreq  = ef * (Ddn_ + (Dup_ - Ddn_) * Math::sigmoid(ef));

        f[VMEAS]      = -vmeas_dot + (vctrl - vmeas) / Tfltr_;
        f[QMEAS]      = -qmeas_dot + (q - qmeas) / Tfltr_;
        f[XQPI]       = -xqpi_dot + sfrz * Math::antiwindup(qpi, Ki_ * erqlim, Qmin_, Qmax_);
        f[XQLAG]      = -xqlag_dot + (qpi - xqlag) / Tfv_;
        f[PMEAS]      = -pmeas_dot + (p - pmeas) / Tp_;
        f[XPPI]       = -xppi_dot + Math::antiwindup(ppi, Kig_ * eplim, Pmin_, Pmax_);
        f[PREF_STATE] = -pref_dot + (ppi - pref) / Tlag_;

        f[V]      = -v * v + vr * vr + vi * vi;
        f[VLDC]   = -vldc * vldc + vldc_r * vldc_r + vldc_i * vldc_i;
        f[VDROOP] = -vdroop + v + Kc_ * q;
        f[VCTRL]  = -vctrl + vcomp_on_ * vldc + vcomp_off_ * vdroop;
        f[SFRZ]   = -sfrz + Math::above(v, Vfrz_);
        f[ERQ]    = -erq + ref_on_ * (vref - vmeas) + ref_off_ * (qref - qmeas);
        f[ERQDB]  = -erqdb + Math::deadband2(erq, dbdlow_, dbdupper_);
        f[ERQLIM] = -erqlim + Math::clamp(erqdb, emin_, emax_);
        f[QPI]    = -qpi + Math::clamp(Kp_ * erqlim + xqpi, Qmin_, Qmax_);
        f[QEXT]   = -Tfv_ * (qext - xqlag) + Tft_ * (qpi - xqlag);

        f[EF]    = -ef + Math::deadband2(freqref - freq, fdbd1_, fdbd2_);
        f[EP]    = -ep + pref_in - pmeas + pfreq;
        f[EPLIM] = -eplim + Math::clamp(ep, femin_, femax_);
        f[PPI]   = -ppi + Math::clamp(Kpg_ * eplim + xppi, Pmin_, Pmax_);
        f[PEXT]  = -pext + freq_on_ * pref;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateInternalResidual()
      {
        y_ext_[4] = ONE<RealT>;
        for (size_t n = 0; n < reference_.size(); ++n)
          y_ext_[5 + n] = static_cast<ScalarT>(reference_[n]);
        this->gatherExternalVariables();
        evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initializeMonitor()
      {
        using Mon = typename ModelDataT::MonitorableVariables;
        monitor_->set(Mon::qext, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::QEXT)]; });
        monitor_->set(Mon::pext, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::PEXT)]; });
        monitor_->set(Mon::vmeas, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::VMEAS)]; });
        monitor_->set(Mon::qmeas, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::QMEAS)]; });
        monitor_->set(Mon::pmeas, [this]
                      { return y_.getData()[static_cast<size_t>(RepcaInternalVariables::PMEAS)]; });
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Repca<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      bool Repca<scalar_type, index_type>::invertClamp(RealT output, RealT lower, RealT upper, RealT& input) const
      {
        const RealT value = static_cast<RealT>(output);

        if (!std::isfinite(value)
            || !std::isfinite(lower)
            || !std::isfinite(upper)
            || lower > upper
            || value < lower
            || value > upper)
        {
          return false;
        }

        const RealT width = upper - lower;
        if (width <= INITIALIZATION_TOLERANCE)
        {
          input = static_cast<RealT>(lower);
          return true;
        }

        const RealT distance_from_lower = value - lower;
        const RealT distance_from_upper = upper - value;
        if (distance_from_lower <= INITIALIZATION_TOLERANCE)
        {
          input = static_cast<RealT>(lower - INITIALIZATION_LIMIT_OFFSET);
          return true;
        }
        if (distance_from_upper <= INITIALIZATION_TOLERANCE)
        {
          input = static_cast<RealT>(upper + INITIALIZATION_LIMIT_OFFSET);
          return true;
        }

        const RealT mu                    = Math::MU<RealT>;
        const RealT scaled_lower_distance = mu * distance_from_lower;
        const RealT scaled_upper_distance = mu * distance_from_upper;
        const RealT log_lower             = logOneMinusExp(scaled_lower_distance);
        const RealT log_upper             = logOneMinusExp(scaled_upper_distance);
        const RealT correction            = (scaled_lower_distance + log_lower - log_upper) / mu;

        input = static_cast<RealT>(lower + correction);
        return std::isfinite(static_cast<RealT>(input));
      }

      template <typename scalar_type, typename index_type>
      bool Repca<scalar_type, index_type>::invertDeadband(RealT output, RealT lower, RealT upper, RealT& input) const
      {
        const RealT value = static_cast<RealT>(output);

        if (!std::isfinite(value)
            || !std::isfinite(lower)
            || !std::isfinite(upper)
            || lower > upper)
        {
          return false;
        }

        const RealT midpoint = HALF<RealT> * lower + HALF<RealT> * upper;
        if (std::abs(value) <= INITIALIZATION_TOLERANCE)
        {
          input = static_cast<RealT>(midpoint);
          return true;
        }

        RealT lower_input = midpoint;
        RealT upper_input = midpoint;
        if (value < ZERO<RealT>)
        {
          lower_input = lower + value;
        }
        else
        {
          upper_input = upper + value;
        }

        const RealT lower_output = Math::deadband2(lower_input, lower, upper);
        const RealT upper_output = Math::deadband2(upper_input, lower, upper);
        if (!std::isfinite(lower_output)
            || !std::isfinite(upper_output)
            || lower_output - value > INITIALIZATION_TOLERANCE
            || value - upper_output > INITIALIZATION_TOLERANCE)
        {
          return false;
        }

        for (std::size_t iteration = 0; iteration < 128; ++iteration)
        {
          const RealT mid = lower_input + HALF<RealT> * (upper_input - lower_input);
          if (Math::deadband2(mid, lower, upper) < value)
          {
            lower_input = mid;
          }
          else
          {
            upper_input = mid;
          }
        }

        const RealT result = lower_input + HALF<RealT> * (upper_input - lower_input);
        input              = static_cast<RealT>(result);
        return std::isfinite(result)
               && std::abs(Math::deadband2(result, lower, upper) - value)
                      <= INITIALIZATION_TOLERANCE;
      }

      template <typename scalar_type, typename index_type>
      typename Repca<scalar_type, index_type>::RealT
      Repca<scalar_type, index_type>::logOneMinusExp(RealT x)
      {
        static constexpr auto log_two = std::numbers::ln2_v<RealT>;

        if (x < log_two)
        {
          return log_two - HALF<RealT> * x
                 + std::log(std::sinh(HALF<RealT> * x));
        }
        return std::log1p(-std::exp(-x));
      }

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

#pragma once

#include <format>
#include <iostream>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/GenrouData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Constructor for a GENROU generator model with saturation
     *
     * Arguments passed to ModelEvaluatorImpl:
     * - Number of equations = 0
     * - Number of independent variables = 0
     * - Number of quadratures = 0
     * - Number of optimization parameters = 0
     */
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(BusT* bus)
      : bus_(bus)
    {
      size_ = 17;
      setDerivedParameters();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(BusT* bus,
                                            RealT p0,
                                            RealT q0,
                                            RealT H,
                                            RealT D,
                                            RealT Ra,
                                            RealT Tdop,
                                            RealT Tdopp,
                                            RealT Tqopp,
                                            RealT Tqop,
                                            RealT Xd,
                                            RealT Xdp,
                                            RealT Xdpp,
                                            RealT Xq,
                                            RealT Xqp,
                                            RealT Xqpp,
                                            RealT Xl,
                                            RealT S10,
                                            RealT S12)
      : bus_(bus),
        p0_(p0),
        q0_(q0),
        H_(H),
        D_(D),
        Ra_(Ra),
        Tdop_(Tdop),
        Tdopp_(Tdopp),
        Tqopp_(Tqopp),
        Tqop_(Tqop),
        Xd_(Xd),
        Xdp_(Xdp),
        Xdpp_(Xdpp),
        Xq_(Xq),
        Xqp_(Xqp),
        Xqpp_(Xqpp),
        Xl_(Xl),
        S10_(S10),
        S12_(S12)
    {
      size_ = 17;
      setDerivedParameters();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(BusT*             bus,
                                            const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();

      size_ = 17;
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(BusT* bus, SignalNodeT* omega, SignalNodeT* pmech, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      ports_.in.template port<GenrouSignalInputs::pmech>().connect(pmech);
      ports_.out.template port<GenrouSignalOutputs::speed>().connect(omega);
      initializeParameters(data);
      initializeMonitor();

      size_ = 17;
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(BusT* bus, SignalNodeT* omega, SignalNodeT* pmech, SignalNodeT* efd, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      ports_.in.template port<GenrouSignalInputs::pmech>().connect(pmech);
      ports_.in.template port<GenrouSignalInputs::efd>().connect(efd);
      ports_.out.template port<GenrouSignalOutputs::speed>().connect(omega);
      initializeParameters(data);
      initializeMonitor();

      size_ = 17;
    }

    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::~Genrou()
    {
    }

    /**
     * @brief Read the parameters out of the model data
     *
     * Every parameter is optional and keeps its documented default when
     * omitted. A non-numeric or non-finite value throws
     * std::invalid_argument. Integer JSON values are accepted for real
     * parameters.
     *
     * @param[in] data Parameters, terminal bus, and monitored-variable selections.
     */
    template <typename scalar_type, typename index_type>
    void Genrou<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = typename ModelDataT::Parameters;
      using Buses  = typename ModelDataT::Buses;

      const Model::ParameterReader reader(data, "Genrou");

      p0_       = reader.readReal(Params::p0).value_or(p0_);
      q0_       = reader.readReal(Params::q0).value_or(q0_);
      H_        = reader.readReal(Params::H).value_or(H_);
      D_        = reader.readReal(Params::D).value_or(D_);
      Ra_       = reader.readReal(Params::Ra).value_or(Ra_);
      Tdop_     = reader.readReal(Params::Tdop).value_or(Tdop_);
      Tdopp_    = reader.readReal(Params::Tdopp).value_or(Tdopp_);
      Tqop_     = reader.readReal(Params::Tqop).value_or(Tqop_);
      Tqopp_    = reader.readReal(Params::Tqopp).value_or(Tqopp_);
      Xd_       = reader.readReal(Params::Xd).value_or(Xd_);
      Xdp_      = reader.readReal(Params::Xdp).value_or(Xdp_);
      Xdpp_     = reader.readReal(Params::Xdpp).value_or(Xdpp_);
      Xq_       = reader.readReal(Params::Xq).value_or(Xq_);
      Xqp_      = reader.readReal(Params::Xqp).value_or(Xqp_);
      Xqpp_     = reader.readReal(Params::Xqpp).value_or(Xqpp_);
      Xl_       = reader.readReal(Params::Xl).value_or(Xl_);
      S10_      = reader.readReal(Params::S10).value_or(S10_);
      S12_      = reader.readReal(Params::S12).value_or(S12_);
      mva_base_ = reader.readReal(Params::mva).value_or(mva_base_);

      if (data.buses.contains(Buses::bus))
      {
        bus_id_ = data.buses.at(Buses::bus);
      }

      setDerivedParameters();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Genrou<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Genrou<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      // Convert monitored terminal values to system base.
      monitor_->set(Variable::ir, [this]
                    { return this->toSystemBase(y_.getData()[15]); });
      monitor_->set(Variable::ii, [this]
                    { return this->toSystemBase(y_.getData()[16]); });
      monitor_->set(Variable::p, [this]
                    { return this->toSystemBase(Vr() * y_.getData()[15] + Vi() * y_.getData()[16]); });
      monitor_->set(Variable::q, [this]
                    { return this->toSystemBase(Vi() * y_.getData()[15] - Vr() * y_.getData()[16]); });
      monitor_->set(Variable::delta, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::omega, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::speed, [this]
                    { return 1.0 + y_.getData()[1]; });
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      auto size = static_cast<size_t>(size_);

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize bus data
      wb_.resize(2);
      h_.resize(2);

      // Resize signal variable data
      ws_.resize(2);
      ws_indices_.resize(2);
      ws_indices_[0] = INVALID_INDEX<IdxT>;
      ws_indices_[1] = INVALID_INDEX<IdxT>;

      // Set output signals
      if (auto speed_port = ports_.out.template port<GenrouSignalOutputs::speed>())
      {
        speed_port.link(&y_.getData()[1], &(this->getVariableIndex(1)));
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Validate the GENROU configuration
     *
     * Checks the terminal bus, the parameters the machine equations divide
     * by, the saturation curve, and attached external signals.
     *
     * @return The configuration checks; passed() when valid.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Genrou<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Genrou");

      checks.check(bus_ != nullptr, "bus pointer is null");
      checks.check(H_ > ZERO<RealT>, "H must be positive");
      checks.check(Tdop_ > ZERO<RealT>, "Tdop must be positive");
      checks.check(Tdopp_ > ZERO<RealT>, "Tdopp must be positive");
      checks.check(Tqop_ > ZERO<RealT>, "Tqop must be positive");
      checks.check(Tqopp_ > ZERO<RealT>, "Tqopp must be positive");
      checks.check(Xd_ != Xl_, "Xd must differ from Xl");
      checks.check(Xdp_ != Xl_, "Xdp must differ from Xl");
      checks.check(Xqp_ != Xl_, "Xqp must differ from Xl");
      checks.check(Ra_ * Ra_ + Xqpp_ * Xqpp_ > ZERO<RealT>, "Ra and Xqpp cannot both be zero");
      checks.check(S10_ >= ZERO<RealT>, "S10 must be non-negative");
      checks.check(S12_ >= ZERO<RealT>, "S12 must be non-negative");
      checks.check(S10_ == ZERO<RealT> || S12_ > ZERO<RealT>,
                   "S12 must be positive when S10 is nonzero");
      checks.check(mva_base_ > ZERO<RealT>, "mva must be positive");

      // An attached port must resolve to readable signal storage. The
      // enumerator is a template argument, so each port names itself once.
      auto check_attached_signal =
          [&]<GenrouSignalInputs input>(const char* name)
      {
        if (ports_.in.template port<input>().connected()
            && !ports_.in.template port<input>().linked())
        {
          checks.fail(std::format("{} signal attached with no linked source", name));
        }
      };

      check_attached_signal.template operator()<GenrouSignalInputs::pmech>("pmech");
      check_attached_signal.template operator()<GenrouSignalInputs::efd>("efd");

      return checks;
    }

    /**
     * Initialization of the generator model
     *
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::initialize()
    {
      if (!this->verifyAndLog())
      {
        Log::error() << "Genrou: cannot initialize with invalid configuration\n";
        return 1;
      }

      // Network Frame Terminal Values
      ScalarT vr  = Vr();
      ScalarT vi  = Vi();
      ScalarT p   = this->toComponentBase(static_cast<ScalarT>(p0_));
      ScalarT q   = this->toComponentBase(static_cast<ScalarT>(q0_));
      ScalarT vm2 = vr * vr + vi * vi;
      ScalarT ir  = (p * vr + q * vi) / vm2;
      ScalarT ii  = (p * vi - q * vr) / vm2;

      // The subtransient-flux magnitude is invariant under the rotor-frame
      // rotation, so saturation is available directly from network quantities.
      const ScalarT Vint_r = vr + Ra_ * ir - Xqpp_ * ii;
      const ScalarT Vint_i = vi + Ra_ * ii + Xqpp_ * ir;
      ScalarT       psipp  = std::sqrt(Vint_r * Vint_r + Vint_i * Vint_i);
      ScalarT       ksat   = SB_ * Math::qramp(psipp - SA_);

      const ScalarT ksat_prime = ONE<RealT> + Xqd_ * ksat;
      const ScalarT xsat_delta = ksat_prime * Xdpp_ + Xq_ - Xqpp_;

      ScalarT delta = std::atan2((vi + Ra_ * ii) * ksat_prime + xsat_delta * ir,
                                 (vr + Ra_ * ir) * ksat_prime - xsat_delta * ii);
      ScalarT id    = ir * std::sin(delta) - ii * std::cos(delta);
      ScalarT iq    = ir * std::cos(delta) + ii * std::sin(delta);
      ScalarT vq    = vr * std::cos(delta) + vi * std::sin(delta) + id * Xqpp_ + iq * Ra_;
      ScalarT Edp   = (Xq1_ - Xqd_ * (Xqp_ - Xqpp_) * ksat) * iq / ksat_prime;
      ScalarT psiqp = Edp + Xq2_ * iq;
      ScalarT psidp = vq - (Xdpp_ - Xl_) * id;
      ScalarT Eqp   = psidp + (Xdp_ - Xl_) * id;

      ScalarT omega(0.0);
      auto*   y  = y_.getData();
      auto*   yp = yp_.getData();

      y[0]           = delta;
      y[1]           = omega;
      y[2]           = Eqp;
      y[3]           = psidp;
      y[4]           = psiqp;
      y[5]           = Edp;
      ScalarT psiqpp = -psiqp * Xq4_ - Edp * Xq5_;
      ScalarT psidpp = psidp * Xd4_ + Eqp * Xd5_;
      y[6]           = psiqpp;
      y[7]           = psidpp;
      y[8] = psipp = std::sqrt(psiqpp * psiqpp + psidpp * psidpp);
      y[9] = ksat = SB_ * Math::qramp(psipp - SA_);
      ScalarT vd  = -psiqpp * (ONE<RealT> + omega);
      vq          = psidpp * (ONE<RealT> + omega);
      y[10]       = vd;
      y[11]       = vq;
      y[12]       = (psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id;
      y[13]       = id;
      y[14]       = iq;
      y[15]       = ir;
      y[16]       = ii;

      // TODO: Writing to an input signal seems very wrong. Check to see if this
      //       could be handled differently

      // Convert Te to system base for governor PM signal.
      ScalarT Te = y[12];
      pmech_set_ = static_cast<RealT>(this->toSystemBase(Te));
      if (auto pmech_port = ports_.in.template port<GenrouSignalInputs::pmech>())
      {
        pmech_port.writeValue(pmech_set_);
      }

      efd_set_ = static_cast<RealT>(Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat);
      if (auto efd_port = ports_.in.template port<GenrouSignalInputs::efd>())
      {
        efd_port.writeValue(efd_set_);
      }

      for (IdxT i = 0; i < size_; ++i)
      {
        yp[static_cast<size_t>(i)] = 0.0;
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      // For DependencyTracking::Variable, set variable numbers
      if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
      {
        this->initializeDependencyTrackingVariableNumbers();
      }

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[static_cast<size_t>(i)] = i < 6;
      }
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Genrou<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y,
        const ScalarT* yp,
        const ScalarT* wb,
        const ScalarT* ws,
        ScalarT*       f)
    {
      /* Read variables */
      ScalarT delta  = y[0];
      ScalarT omega  = y[1];
      ScalarT Eqp    = y[2];
      ScalarT psidp  = y[3];
      ScalarT psiqp  = y[4];
      ScalarT Edp    = y[5];
      ScalarT psiqpp = y[6];
      ScalarT psidpp = y[7];
      ScalarT psipp  = y[8];
      ScalarT ksat   = y[9];
      ScalarT vd     = y[10];
      ScalarT vq     = y[11];
      ScalarT telec  = y[12];
      ScalarT id     = y[13];
      ScalarT iq     = y[14];
      ScalarT ir     = y[15];
      ScalarT ii     = y[16];

      /* Read derivatives */
      ScalarT delta_dot = yp[0];
      ScalarT omega_dot = yp[1];
      ScalarT Eqp_dot   = yp[2];
      ScalarT psidp_dot = yp[3];
      ScalarT psiqp_dot = yp[4];
      ScalarT Edp_dot   = yp[5];

      // Set coupling variable aliases
      ScalarT vr = wb[0];
      ScalarT vi = wb[1];

      // Set signal variable aliases
      ScalarT pmech = this->toComponentBase(ws[0]);
      ScalarT efd   = ws[1];

      static constexpr auto pi = std::numbers::pi_v<RealT>;

      /* 6 Genrou differential equations */
      f[0] = delta_dot - omega * (TWO<RealT> * pi * freq_system_base_);
      f[1] = omega_dot - (ONE<RealT> / (TWO<RealT> * H_)) * ((pmech - D_ * omega) / (ONE<RealT> + omega) - telec);
      f[2] = Eqp_dot - (ONE<RealT> / Tdop_) * (efd - (Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat));
      f[3] = psidp_dot - (ONE<RealT> / Tdopp_) * (Eqp - psidp - Xd2_ * id);
      f[4] = psiqp_dot - (ONE<RealT> / Tqopp_) * (Edp - psiqp + Xq2_ * iq);
      f[5] = Edp_dot - (ONE<RealT> / Tqop_) * (-Edp + Xqd_ * psiqpp * ksat + Xq1_ * (iq - Xq3_ * (Edp + iq * Xq2_ - psiqp)));

      /* 11 Genrou algebraic equations */
      f[6]  = psiqpp - (-psiqp * Xq4_ - Edp * Xq5_);
      f[7]  = psidpp - (psidp * Xd4_ + Eqp * Xd5_);
      f[8]  = psipp - std::sqrt((psidpp * psidpp) + (psiqpp * psiqpp));
      f[9]  = ksat - SB_ * Math::qramp(psipp - SA_);
      f[10] = vd + psiqpp * (ONE<RealT> + omega);
      f[11] = vq - psidpp * (ONE<RealT> + omega);
      f[12] = telec - ((psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id);
      f[13] = id - (ir * std::sin(delta) - ii * std::cos(delta));
      f[14] = iq - (ir * std::cos(delta) + ii * std::sin(delta));
      f[15] = ir - (G_ * (std::sin(delta) * vd + std::cos(delta) * vq - vr) - B_ * (-std::cos(delta) * vd + std::sin(delta) * vq - vi));
      f[16] = ii - (B_ * (std::sin(delta) * vd + std::cos(delta) * vq - vr) + G_ * (-std::cos(delta) * vd + std::sin(delta) * vq - vi));

      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Genrou<scalar_type, index_type>::evaluateBusResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* wb,
        ScalarT*                        h)
    {
      ScalarT ir = y[15];
      ScalarT ii = y[16];

      // Convert current injection to system base for the network.
      h[0] = this->toSystemBase(ir);
      h[1] = this->toSystemBase(ii);

      return 0;
    }

    /**
     * \brief Residual evaluation and contribution to the connected bus
     *
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::evaluateResidual()
    {
      auto* ws = ws_.getData();

      // Mechanical Power
      ws[0] = pmech_set_;
      if (auto pmech_port = ports_.in.template port<GenrouSignalInputs::pmech>())
      {
        ws[0]          = pmech_port.readSignal();
        ws_indices_[0] = pmech_port.signalVariableIndex();
      }

      // Exciter Efield
      ws[1] = efd_set_;
      if (auto efd_port = ports_.in.template port<GenrouSignalInputs::efd>())
      {
        ws[1]          = efd_port.readSignal();
        ws_indices_[1] = efd_port.signalVariableIndex();
      }

      // Bus voltages
      auto* wb = wb_.getData();
      wb[0]    = Vr();
      wb[1]    = Vi();

      // Residual evaluation
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      auto*       h  = h_.getData();
      evaluateInternalResidual(y, yp, wb, ws, f);
      evaluateBusResidual(y, yp, wb, h);

      // Genrou contribution to bus algebraic equations
      Ir() += h[0];
      Ii() += h[1];

      if (bus_->size() > 0)
      {
        bus_->getResidual().setDataUpdated();
      }
      f_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Genrou<scalar_type, index_type>::setDerivedParameters()
    {
      SA_ = 0;
      SB_ = 0;
      if (S12_ != 0)
      {
        RealT s112 = std::sqrt(S10_ / S12_);

        SA_ = (1.2 * s112 + ONE<RealT>) / (s112 + ONE<RealT>);
        SB_ = (1.2 * s112 - ONE<RealT>) / (s112 - ONE<RealT>);
        if (SB_ < SA_)
        {
          SA_ = SB_;
        }
        SB_ = S12_ / ((SA_ - 1.2) * (SA_ - 1.2));
      }
      Xd1_ = Xd_ - Xdp_;
      Xd2_ = Xdp_ - Xl_;
      Xd3_ = (Xdp_ - Xdpp_) / (Xd2_ * Xd2_);
      Xd4_ = (Xdp_ - Xdpp_) / Xd2_;
      Xd5_ = (Xdpp_ - Xl_) / Xd2_;
      Xq1_ = Xq_ - Xqp_;
      Xq2_ = Xqp_ - Xl_;
      Xq3_ = (Xqp_ - Xqpp_) / (Xq2_ * Xq2_);
      Xq4_ = (Xqp_ - Xqpp_) / Xq2_;
      Xq5_ = (Xqpp_ - Xl_) / Xq2_;
      Xqd_ = (Xq_ - Xl_) / (Xd_ - Xl_);
      G_   = Ra_ / (Ra_ * Ra_ + Xqpp_ * Xqpp_);
      B_   = -Xqpp_ / (Ra_ * Ra_ + Xqpp_ * Xqpp_);
      this->setComponentBase(mva_base_ * static_cast<RealT>(1.0e6));
    }
  } // namespace PhasorDynamics
} // namespace GridKit

#pragma once

#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
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
      : bus_(bus),
        bus_id_(0),
        p0_(0.),
        q0_(0.),
        H_(3.),
        D_(0.),
        Ra_(0.),
        Tdop_(7.),
        Tdopp_(.04),
        Tqopp_(.05),
        Tqop_(.75),
        Xd_(2.1),
        Xdp_(0.2),
        Xdpp_(0.18),
        Xq_(.5),
        Xqp_(.5),
        Xqpp_(.18),
        Xl_(.15),
        S10_(0.),
        S12_(0.),
        mva_base_(100.)
    {
      size_ = 15;
      setDerivedParams();
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
        bus_id_(0),
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
        S12_(S12),
        mva_base_(100.)
    {
      size_ = 15;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(BusT* bus, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();

      size_ = 15;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(BusT* bus, SignalT* omega, SignalT* pmech, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      signals_.template attachSignalNode<GenrouExternalVariables::PM>(pmech);
      signals_.template assignSignalNode<GenrouInternalVariables::OMEGA>(omega);
      initializeParameters(data);
      initializeMonitor();

      size_ = 15;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(BusT* bus, SignalT* omega, SignalT* pmech, SignalT* efd, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      signals_.template attachSignalNode<GenrouExternalVariables::PM>(pmech);
      signals_.template assignSignalNode<GenrouInternalVariables::OMEGA>(omega);
      signals_.template attachSignalNode<GenrouExternalVariables::EFD>(efd);
      initializeParameters(data);
      initializeMonitor();

      size_ = 15;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::~Genrou()
    {
    }

    /// Helper function to extract and assign model parameters from the model's associated
    /// data structure.
    template <typename scalar_type, typename index_type>
    void Genrou<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      using Buses     = typename ModelDataT::Buses;
      if (data.parameters.contains(Parameter::p0))
      {
        p0_ = std::get<RealT>(data.parameters.at(Parameter::p0));
      }

      if (data.parameters.contains(Parameter::q0))
      {
        q0_ = std::get<RealT>(data.parameters.at(Parameter::q0));
      }

      if (data.parameters.contains(Parameter::H))
      {
        H_ = std::get<RealT>(data.parameters.at(Parameter::H));
      }

      if (data.parameters.contains(Parameter::D))
      {
        D_ = std::get<RealT>(data.parameters.at(Parameter::D));
      }

      if (data.parameters.contains(Parameter::Ra))
      {
        Ra_ = std::get<RealT>(data.parameters.at(Parameter::Ra));
      }

      if (data.parameters.contains(Parameter::Tdop))
      {
        Tdop_ = std::get<RealT>(data.parameters.at(Parameter::Tdop));
      }

      if (data.parameters.contains(Parameter::Tdopp))
      {
        Tdopp_ = std::get<RealT>(data.parameters.at(Parameter::Tdopp));
      }

      if (data.parameters.contains(Parameter::Tqopp))
      {
        Tqopp_ = std::get<RealT>(data.parameters.at(Parameter::Tqopp));
      }

      if (data.parameters.contains(Parameter::Tqop))
      {
        Tqop_ = std::get<RealT>(data.parameters.at(Parameter::Tqop));
      }

      if (data.parameters.contains(Parameter::Xd))
      {
        Xd_ = std::get<RealT>(data.parameters.at(Parameter::Xd));
      }

      if (data.parameters.contains(Parameter::Xdp))
      {
        Xdp_ = std::get<RealT>(data.parameters.at(Parameter::Xdp));
      }

      if (data.parameters.contains(Parameter::Xdpp))
      {
        Xdpp_ = std::get<RealT>(data.parameters.at(Parameter::Xdpp));
      }

      if (data.parameters.contains(Parameter::Xq))
      {
        Xq_ = std::get<RealT>(data.parameters.at(Parameter::Xq));
      }

      if (data.parameters.contains(Parameter::Xqp))
      {
        Xqp_ = std::get<RealT>(data.parameters.at(Parameter::Xqp));
      }

      if (data.parameters.contains(Parameter::Xqpp))
      {
        Xqpp_ = std::get<RealT>(data.parameters.at(Parameter::Xqpp));
      }

      if (data.parameters.contains(Parameter::Xl))
      {
        Xl_ = std::get<RealT>(data.parameters.at(Parameter::Xl));
      }

      if (data.parameters.contains(Parameter::S10))
      {
        S10_ = std::get<RealT>(data.parameters.at(Parameter::S10));
      }

      if (data.parameters.contains(Parameter::S12))
      {
        S12_ = std::get<RealT>(data.parameters.at(Parameter::S12));
      }

      if (data.parameters.contains(Parameter::mva))
      {
        mva_base_ = std::get<RealT>(data.parameters.at(Parameter::mva));
      }

      if (data.buses.contains(Buses::bus))
      {
        bus_id_ = data.buses.at(Buses::bus);
      }
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Genrou<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    // System base -> machine base when reading system values.
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::ScalarT Genrou<scalar_type, index_type>::toMachineBase(ScalarT value) const
    {
      return value * va_system_base_ / va_machine_base_;
    }

    // Machine base -> system base for network and signal output.
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::ScalarT Genrou<scalar_type, index_type>::toSystemBase(ScalarT value) const
    {
      return value / toMachineBase(static_cast<ScalarT>(ONE<RealT>));
    }

    template <typename scalar_type, typename index_type>
    void Genrou<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      // Convert monitored terminal values to system base.
      monitor_->set(Variable::ir, [this]
                    { return toSystemBase(y_.getData()[13]); });
      monitor_->set(Variable::ii, [this]
                    { return toSystemBase(y_.getData()[14]); });
      monitor_->set(Variable::p, [this]
                    { return toSystemBase(Vr() * y_.getData()[13] + Vi() * y_.getData()[14]); });
      monitor_->set(Variable::q, [this]
                    { return toSystemBase(Vi() * y_.getData()[13] - Vr() * y_.getData()[14]); });
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
      if (signals_.template isAssigned<GenrouInternalVariables::OMEGA>())
      {
        auto* y = y_.getData();
        signals_.template getSignalNode<GenrouInternalVariables::OMEGA>()->set(&y[1], &(this->getVariableIndex(1)));
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::verify() const
    {
      static constexpr auto PM  = GenrouExternalVariables::PM;
      static constexpr auto EFD = GenrouExternalVariables::EFD;

      int ret = 0;

      if (signals_.template isAttached<PM>())
      {
        if (!signals_.template isLinked<PM>())
        {
          Log::error() << "Genrou: pmech signal attached with no linked governor\n";
          ret += 1;
        }
      }

      if (signals_.template isAttached<EFD>())
      {
        if (!signals_.template isLinked<EFD>())
        {
          Log::error() << "Genrou: efd signal attached with no linked exciter\n";
          ret += 1;
        }
      }

      return ret;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::initialize()
    {
      // Network Frame Terminal Values
      ScalarT vr  = Vr();
      ScalarT vi  = Vi();
      ScalarT p   = toMachineBase(static_cast<ScalarT>(p0_));
      ScalarT q   = toMachineBase(static_cast<ScalarT>(q0_));
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
      y[10]       = (psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id;
      y[11]       = id;
      y[12]       = iq;
      y[13]       = ir;
      y[14]       = ii;

      ScalarT Te = y[10];
      // Convert Te to system base for governor PM signal.
      pmech_set_ = toSystemBase(Te);
      if (signals_.template isAttached<GenrouExternalVariables::PM>())
      {
        signals_.template writeExternalVariable<GenrouExternalVariables::PM>(pmech_set_);
      }

      efd_set_ = Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat;
      if (signals_.template isAttached<GenrouExternalVariables::EFD>())
      {
        signals_.template writeExternalVariable<GenrouExternalVariables::EFD>(efd_set_);
      }

      for (IdxT i = 0; i < size_; ++i)
      {
        yp[static_cast<size_t>(i)] = 0.0;
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

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
      ScalarT telec  = y[10];
      ScalarT id     = y[11];
      ScalarT iq     = y[12];
      ScalarT ir     = y[13];
      ScalarT ii     = y[14];

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
      ScalarT pmech = toMachineBase(ws[0]);
      ScalarT efd   = ws[1];

      static constexpr auto pi = std::numbers::pi_v<RealT>;

      // Set Rotor Angle commputation
      const ScalarT sin_delta = std::sin(delta);
      const ScalarT cos_delta = std::cos(delta);

      // Internal Voltage
      const ScalarT Vint_r = (-sin_delta * psiqpp + cos_delta * psidpp) * (ONE<RealT> + omega);
      const ScalarT Vint_i = (cos_delta * psiqpp + sin_delta * psidpp) * (ONE<RealT> + omega);

      /* 6 Genrou differential equations */
      f[0] = delta_dot - omega * (TWO<RealT> * pi * freq_system_base_);
      f[1] = omega_dot - (ONE<RealT> / (TWO<RealT> * H_)) * ((pmech - D_ * omega) / (ONE<RealT> + omega) - telec);
      f[2] = Eqp_dot - (ONE<RealT> / Tdop_) * (efd - (Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat));
      f[3] = psidp_dot - (ONE<RealT> / Tdopp_) * (Eqp - psidp - Xd2_ * id);
      f[4] = psiqp_dot - (ONE<RealT> / Tqopp_) * (Edp - psiqp + Xq2_ * iq);
      f[5] = Edp_dot - (ONE<RealT> / Tqop_) * (-Edp + Xqd_ * psiqpp * ksat + Xq1_ * (iq - Xq3_ * (Edp + iq * Xq2_ - psiqp)));

      /* 7 Genrou algebraic equations */
      f[6]  = psiqpp - (-psiqp * Xq4_ - Edp * Xq5_);
      f[7]  = psidpp - (psidp * Xd4_ + Eqp * Xd5_);
      f[8]  = psipp - std::sqrt((psidpp * psidpp) + (psiqpp * psiqpp));
      f[9]  = ksat - SB_ * Math::qramp(psipp - SA_);
      f[10] = telec - ((psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id);
      f[11] = id - (ir * sin_delta - ii * cos_delta);
      f[12] = iq - (ir * cos_delta + ii * sin_delta);

      /* 2 Genrou network equations */
      f[13] = ir + G_ * vr - B_ * vi - (G_ * Vint_r - B_ * Vint_i);
      f[14] = ii + B_ * vr + G_ * vi - (B_ * Vint_r + G_ * Vint_i);

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
      ScalarT ir = y[13];
      ScalarT ii = y[14];

      // Convert current injection to system base for the network.
      h[0] = toSystemBase(ir);
      h[1] = toSystemBase(ii);

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
      if (signals_.template isAttached<GenrouExternalVariables::PM>())
      {
        ws[0]          = signals_.template readExternalVariable<GenrouExternalVariables::PM>();
        ws_indices_[0] = signals_.template readExternalVariableIndex<GenrouExternalVariables::PM>();
      }

      // Exciter Efield
      ws[1] = efd_set_;
      if (signals_.template isAttached<GenrouExternalVariables::EFD>())
      {
        ws[1]          = signals_.template readExternalVariable<GenrouExternalVariables::EFD>();
        ws_indices_[1] = signals_.template readExternalVariableIndex<GenrouExternalVariables::EFD>();
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

      f_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Genrou<scalar_type, index_type>::setDerivedParams()
    {
      SA_ = 0;
      SB_ = 0;
      if (S12_ != 0)
      {
        RealT s112 = std::sqrt(S10_ / S12_);

        SA_ = (1.2 * s112 + ONE<RealT>) / (s112 + ONE<RealT>);
        SB_ = (1.2 * s112 - ONE<RealT>) / (s112 - ONE<RealT>);
        if (SB_ < SA_)
          SA_ = SB_;
        SB_ = S12_ / ((SA_ - 1.2) * (SA_ - 1.2));
      }
      Xd1_             = Xd_ - Xdp_;
      Xd2_             = Xdp_ - Xl_;
      Xd3_             = (Xdp_ - Xdpp_) / (Xd2_ * Xd2_);
      Xd4_             = (Xdp_ - Xdpp_) / Xd2_;
      Xd5_             = (Xdpp_ - Xl_) / Xd2_;
      Xq1_             = Xq_ - Xqp_;
      Xq2_             = Xqp_ - Xl_;
      Xq3_             = (Xqp_ - Xqpp_) / (Xq2_ * Xq2_);
      Xq4_             = (Xqp_ - Xqpp_) / Xq2_;
      Xq5_             = (Xqpp_ - Xl_) / Xq2_;
      Xqd_             = (Xq_ - Xl_) / (Xd_ - Xl_);
      G_               = Ra_ / (Ra_ * Ra_ + Xqpp_ * Xqpp_);
      B_               = -Xqpp_ / (Ra_ * Ra_ + Xqpp_ * Xqpp_);
      va_machine_base_ = mva_base_ * static_cast<RealT>(1.0e6);
    }
  } // namespace PhasorDynamics
} // namespace GridKit

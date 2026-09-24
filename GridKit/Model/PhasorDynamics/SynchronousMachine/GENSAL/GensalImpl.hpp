#pragma once

#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/Gensal.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/GensalData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Constructor for a GENSAL generator model with saturation
     */
    template <typename scalar_type, typename index_type>
    Gensal<scalar_type, index_type>::Gensal(BusT* bus, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();

      size_ = 5;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    Gensal<scalar_type, index_type>::~Gensal()
    {
    }

    /// Helper function to extract and assign model parameters from the model's associated
    /// data structure.
    template <typename scalar_type, typename index_type>
    void Gensal<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
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
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Gensal<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Gensal<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      // Terminal quantities are evaluated from the states, not read back from
      // them, and are converted to system base for reporting.
      monitor_->set(Variable::ir, [this]
                    { return this->toSystemBase(algebraicState().ir); });
      monitor_->set(Variable::ii, [this]
                    { return this->toSystemBase(algebraicState().ii); });
      monitor_->set(Variable::p,
                    [this]
                    {
                      const AlgebraicState s = algebraicState();
                      return this->toSystemBase(Vr() * s.ir + Vi() * s.ii);
                    });
      monitor_->set(Variable::q,
                    [this]
                    {
                      const AlgebraicState s = algebraicState();
                      return this->toSystemBase(Vi() * s.ir - Vr() * s.ii);
                    });
      monitor_->set(Variable::delta, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::omega, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::speed, [this]
                    { return 1.0 + y_.getData()[1]; });
      monitor_->set(Variable::Eqp, [this]
                    { return y_.getData()[2]; });
      monitor_->set(Variable::psidp, [this]
                    { return y_.getData()[3]; });
      monitor_->set(Variable::psiqpp, [this]
                    { return y_.getData()[4]; });
      monitor_->set(Variable::psidpp, [this]
                    { return algebraicState().psidpp; });
      monitor_->set(Variable::vd, [this]
                    { return algebraicState().vd; });
      monitor_->set(Variable::vq, [this]
                    { return algebraicState().vq; });
      monitor_->set(Variable::te, [this]
                    { return algebraicState().telec; });
      monitor_->set(Variable::id, [this]
                    { return algebraicState().id; });
      monitor_->set(Variable::iq, [this]
                    { return algebraicState().iq; });
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::allocate()
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
      if (auto speed_port = ports_.out.template port<GensalSignalOutputs::speed>())
      {
        speed_port.link(&y_.getData()[1], &(this->getVariableIndex(1)));
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::verify() const
    {
      int ret = 0;

      auto pmech_port = ports_.in.template port<GensalSignalInputs::pmech>();
      if (pmech_port.connected() && !pmech_port.linked())
      {
        Log::error() << "Gensal: pmech signal attached with no linked governor\n";
        ret += 1;
      }

      auto efd_port = ports_.in.template port<GensalSignalInputs::efd>();
      if (efd_port.connected() && !efd_port.linked())
      {
        Log::error() << "Gensal: efd signal attached with no linked exciter\n";
        ret += 1;
      }

      return ret;
    }

    /**
     * Initialization of the generator model
     *
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::initialize()
    {
      // Network frame terminal values
      ScalarT vr  = Vr();
      ScalarT vi  = Vi();
      ScalarT p   = this->toComponentBase(static_cast<ScalarT>(p0_));
      ScalarT q   = this->toComponentBase(static_cast<ScalarT>(q0_));
      ScalarT vm2 = vr * vr + vi * vi;
      ScalarT ir  = (p * vr + q * vi) / vm2;
      ScalarT ii  = (p * vi - q * vr) / vm2;

      ScalarT Er    = vr + Ra_ * ir - Xq_ * ii;
      ScalarT Ei    = vi + Ra_ * ii + Xq_ * ir;
      ScalarT delta = std::atan2(Ei, Er);
      ScalarT omega(0.0);

      ScalarT id     = ir * std::sin(delta) - ii * std::cos(delta);
      ScalarT iq     = ir * std::cos(delta) + ii * std::sin(delta);
      ScalarT psiqpp = -Xq2_ * iq;
      ScalarT vq     = vr * std::cos(delta) + vi * std::sin(delta) + id * Xdpp_ + iq * Ra_;
      ScalarT psidpp = vq / (ONE<RealT> + omega);
      ScalarT psidp  = psidpp - (Xdpp_ - Xl_) * id;
      ScalarT Eqp    = psidp + Xd2_ * id;

      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0] = delta;
      y[1] = omega;
      y[2] = Eqp;
      y[3] = psidp;
      y[4] = psiqpp;

      // Algebraic quantities consistent with the converged states use the
      // runtime chain so initial setpoints match residual evaluation exactly.
      const ScalarT        wb[2] = {vr, vi};
      const AlgebraicState s     = evaluateAlgebraicState(y, wb);

      // Convert Te to system base for governor PM signal.
      pmech_set_ = static_cast<RealT>(this->toSystemBase(s.telec));
      if (auto pmech_port = ports_.in.template port<GensalSignalInputs::pmech>())
      {
        pmech_port.writeValue(pmech_set_);
      }

      efd_set_ = static_cast<RealT>(Eqp + Xd1_ * (s.id + Xd3_ * (Eqp - psidp - Xd2_ * s.id)) + Eqp * s.ksat);
      if (auto efd_port = ports_.in.template port<GensalSignalInputs::efd>())
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
     *
     * Every unknown the machine carries is a differential state; its algebraic
     * quantities are evaluated inline rather than solved for.
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::tagDifferentiable()
    {
      tag_.assign(static_cast<size_t>(size_), true);
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
    int Gensal<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Evaluate the machine's algebraic quantities.
     *
     * Subtransient flux, saturation, internal voltage, terminal current,
     * rotor-frame current and electrical torque are an explicit feed-forward
     * chain over the five states and the terminal voltage, with no algebraic
     * loop anywhere in it. They are therefore evaluated here rather than
     * carried as unknowns and solved for.
     *
     * Both residuals and the variable monitor go through this function, so the
     * chain has exactly one definition for Enzyme to differentiate.
     *
     * @param[in] y  - Internal variables
     * @param[in] wb - Bus variables
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline typename Gensal<scalar_type, index_type>::AlgebraicState
    Gensal<scalar_type, index_type>::evaluateAlgebraicState(const ScalarT* y, const ScalarT* wb) const
    {
      /* Read variables */
      const ScalarT delta  = y[0];
      const ScalarT omega  = y[1];
      const ScalarT Eqp    = y[2];
      const ScalarT psidp  = y[3];
      const ScalarT psiqpp = y[4];

      // Set coupling variable aliases
      const ScalarT vr = wb[0];
      const ScalarT vi = wb[1];

      // Set Rotor Angle computation
      const ScalarT sin_delta = std::sin(delta);
      const ScalarT cos_delta = std::cos(delta);

      AlgebraicState s;

      // Subtransient flux linkage and saturation on the transient voltage
      s.psidpp = psidp * Xd4_ + Eqp * Xd5_;
      s.ksat   = SB_ * Math::qramp(Eqp - SA_);

      // Internal voltage in the rotor frame
      s.vd = -psiqpp * (ONE<RealT> + omega);
      s.vq = s.psidpp * (ONE<RealT> + omega);

      // Norton current injection at the terminal
      const ScalarT Vint_r = sin_delta * s.vd + cos_delta * s.vq;
      const ScalarT Vint_i = -cos_delta * s.vd + sin_delta * s.vq;
      s.ir                 = G_ * (Vint_r - vr) - B_ * (Vint_i - vi);
      s.ii                 = B_ * (Vint_r - vr) + G_ * (Vint_i - vi);

      // Rotor-frame currents and the electrical torque they develop
      s.id    = s.ir * sin_delta - s.ii * cos_delta;
      s.iq    = s.ir * cos_delta + s.ii * sin_delta;
      s.telec = (s.psidpp - s.id * Xdpp_) * s.iq - (psiqpp - s.iq * Xdpp_) * s.id;

      return s;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Gensal<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y,
        const ScalarT* yp,
        const ScalarT* wb,
        const ScalarT* ws,
        ScalarT*       f)
    {
      /* Read variables */
      ScalarT omega  = y[1];
      ScalarT Eqp    = y[2];
      ScalarT psidp  = y[3];
      ScalarT psiqpp = y[4];

      /* Read derivatives */
      ScalarT delta_dot  = yp[0];
      ScalarT omega_dot  = yp[1];
      ScalarT Eqp_dot    = yp[2];
      ScalarT psidp_dot  = yp[3];
      ScalarT psiqpp_dot = yp[4];

      // Set signal variable aliases
      ScalarT pmech = this->toComponentBase(ws[0]);
      ScalarT efd   = ws[1];

      static constexpr auto pi = std::numbers::pi_v<RealT>;

      // Algebraic quantities, evaluated rather than solved for
      const AlgebraicState s = evaluateAlgebraicState(y, wb);

      /* 5 Gensal differential equations */
      f[0] = delta_dot - omega * (TWO<RealT> * pi * freq_system_base_);
      f[1] = omega_dot - (ONE<RealT> / (TWO<RealT> * H_)) * ((pmech - D_ * omega) / (ONE<RealT> + omega) - s.telec);
      f[2] = Eqp_dot - (ONE<RealT> / Tdop_) * (efd - (Eqp + Xd1_ * (s.id + Xd3_ * (Eqp - psidp - Xd2_ * s.id)) + Eqp * s.ksat));
      f[3] = psidp_dot - (ONE<RealT> / Tdopp_) * (Eqp - psidp - Xd2_ * s.id);
      f[4] = psiqpp_dot - (ONE<RealT> / Tqopp_) * (-psiqpp - Xq2_ * s.iq);

      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Gensal<scalar_type, index_type>::evaluateBusResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      const AlgebraicState s = evaluateAlgebraicState(y, wb);

      // Convert current injection to system base for the network.
      h[0] = this->toSystemBase(s.ir);
      h[1] = this->toSystemBase(s.ii);

      return 0;
    }

    /**
     * \brief Residual evaluation and contribution to the connected bus
     *
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::evaluateResidual()
    {
      auto* ws = ws_.getData();

      // Mechanical Power
      ws[0] = pmech_set_;
      if (auto pmech_port = ports_.in.template port<GensalSignalInputs::pmech>())
      {
        ws[0]          = pmech_port.readSignal();
        ws_indices_[0] = pmech_port.signalVariableIndex();
      }

      // Exciter Efield
      ws[1] = efd_set_;
      if (auto efd_port = ports_.in.template port<GensalSignalInputs::efd>())
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

      // Gensal contribution to bus algebraic equations
      Ir() += h[0];
      Ii() += h[1];

      f_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Gensal<scalar_type, index_type>::setDerivedParams()
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
      Xq2_ = Xq_ - Xdpp_;
      G_   = Ra_ / (Ra_ * Ra_ + Xdpp_ * Xdpp_);
      B_   = -Xdpp_ / (Ra_ * Ra_ + Xdpp_ * Xdpp_);
      this->setComponentBase(mva_base_ * static_cast<RealT>(1.0e6));
    }
  } // namespace PhasorDynamics
} // namespace GridKit

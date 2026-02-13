#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp> // <- TODO: Temporary, to be removed.
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>
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
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, IdxT unit_id)
      : bus_(bus),
        bus_id_(0),
        unit_id_(unit_id),
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
      size_ = 19;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus,
                                  IdxT      unit_id,
                                  RealT     p0,
                                  RealT     q0,
                                  RealT     H,
                                  RealT     D,
                                  RealT     Ra,
                                  RealT     Tdop,
                                  RealT     Tdopp,
                                  RealT     Tqopp,
                                  RealT     Tqop,
                                  RealT     Xd,
                                  RealT     Xdp,
                                  RealT     Xdpp,
                                  RealT     Xq,
                                  RealT     Xqp,
                                  RealT     Xqpp,
                                  RealT     Xl,
                                  RealT     S10,
                                  RealT     S12)
      : bus_(bus),
        bus_id_(0),
        unit_id_(unit_id),
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
      size_ = 19;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, const model_data_type& data)
      : bus_(bus),
        unit_id_(1),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();

      size_ = 19;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, signal_type* omega, signal_type* pmech, const model_data_type& data)
      : bus_(bus),
        unit_id_(1),
        monitor_(std::make_unique<MonitorT>(data))
    {
      signals_.template attachSignalNode<GenrouExternalVariables::PM>(pmech);
      signals_.template assignSignalNode<GenrouInternalVariables::OMEGA>(omega);
      initializeParameters(data);
      initializeMonitor();

      size_ = 19;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, signal_type* omega, signal_type* pmech, signal_type* efd, const model_data_type& data)
      : bus_(bus),
        unit_id_(1),
        monitor_(std::make_unique<MonitorT>(data))
    {
      signals_.template attachSignalNode<GenrouExternalVariables::PM>(pmech);
      signals_.template assignSignalNode<GenrouInternalVariables::OMEGA>(omega);
      signals_.template attachSignalNode<GenrouExternalVariables::EFD>(efd);
      initializeParameters(data);
      initializeMonitor();

      size_ = 19;
      setDerivedParams();
    }

    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::~Genrou()
    {
    }

    /// Helper function to extract and assign model parameters from the model's associated
    /// data structure.
    template <class ScalarT, typename IdxT>
    void Genrou<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
    {
      if (data.parameters.contains(model_data_type::Parameters::p0))
      {
        p0_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::p0));
      }

      if (data.parameters.contains(model_data_type::Parameters::q0))
      {
        q0_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::q0));
      }

      if (data.parameters.contains(model_data_type::Parameters::H))
      {
        H_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::H));
      }

      if (data.parameters.contains(model_data_type::Parameters::D))
      {
        D_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::D));
      }

      if (data.parameters.contains(model_data_type::Parameters::Ra))
      {
        Ra_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Ra));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tdop))
      {
        Tdop_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tdop));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tdopp))
      {
        Tdopp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tdopp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tqopp))
      {
        Tqopp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tqopp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tqop))
      {
        Tqop_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tqop));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xd))
      {
        Xd_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xd));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xdp))
      {
        Xdp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xdp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xdpp))
      {
        Xdpp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xdpp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xq))
      {
        Xq_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xq));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xqp))
      {
        Xqp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xqp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xqpp))
      {
        Xqpp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xqpp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xl))
      {
        Xl_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xl));
      }

      if (data.parameters.contains(model_data_type::Parameters::S10))
      {
        S10_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::S10));
      }

      if (data.parameters.contains(model_data_type::Parameters::S12))
      {
        S12_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::S12));
      }

      if (data.parameters.contains(model_data_type::Parameters::mva_base))
      {
        mva_base_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::mva_base));
      }

      if (data.ports.contains(model_data_type::Ports::bus))
      {
        bus_id_ = data.ports.at(model_data_type::Ports::bus);
      }
    }

    template <class ScalarT, typename IdxT>
    const Model::VariableMonitorBase* Genrou<ScalarT, IdxT>::getMonitor() const
    {
      return monitor_.get();
    }

    template <class ScalarT, typename IdxT>
    void Genrou<ScalarT, IdxT>::initializeMonitor()
    {
      using Variable = typename model_data_type::MonitorableVariables;
      monitor_->set(Variable::ir, [this]
                    { return y_[15]; });
      monitor_->set(Variable::ii, [this]
                    { return y_[16]; });
      // monitor_->set(Variable::p, [this] { return ?(); });
      // monitor_->set(Variable::q, [this] { return ?(); });
      monitor_->set(Variable::delta, [this]
                    { return y_[0]; });
      monitor_->set(Variable::omega, [this]
                    { return y_[1]; });
    }

    /**
     * @brief Set the component ID
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::allocate()
    {
      // Resize component model data
      auto size = static_cast<size_t>(size_);
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);
      abs_tol_.resize(static_cast<size_t>(size_));

      // Resize bus data
      wb_.resize(2);
      h_.resize(2);

      // Resize signal variable data
      ws_.resize(2);
      ws_indices_.resize(2);
      ws_indices_[0] = INVALID_INDEX<IdxT>;
      ws_indices_[1] = INVALID_INDEX<IdxT>;

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Set output signals
      if (signals_.template isAssigned<GenrouInternalVariables::OMEGA>())
      {
        signals_.template getSignalNode<GenrouInternalVariables::OMEGA>()->set(&y_[1], &(this->getVariableIndex(1)));
      }

      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::verify() const
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
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::initialize()
    {
      /* Initialization tricks -- assuming NO saturation */
      ScalarT vr    = Vr();
      ScalarT vi    = Vi();
      ScalarT p     = static_cast<ScalarT>(p0_) * mva_system_base_ / mva_base_;
      ScalarT q     = static_cast<ScalarT>(q0_) * mva_system_base_ / mva_base_;
      ScalarT vm2   = vr * vr + vi * vi;
      ScalarT Er    = vr + (Ra_ * p * vr + Ra_ * q * vi - Xq_ * p * vi + Xq_ * q * vr) / vm2;
      ScalarT Ei    = vi + (Ra_ * p * vi - Ra_ * q * vr + Xq_ * p * vr + Xq_ * q * vi) / vm2;
      ScalarT delta = std::atan2(Ei, Er);
      ScalarT omega(0.0);
      ScalarT ir     = (p * vr + q * vi) / vm2;
      ScalarT ii     = (p * vi - q * vr) / vm2;
      ScalarT id     = ir * std::sin(delta) - ii * std::cos(delta);
      ScalarT iq     = ir * std::cos(delta) + ii * std::sin(delta);
      ScalarT vd     = vr * std::sin(delta) - vi * std::cos(delta) + id * Ra_ - iq * Xqpp_;
      ScalarT vq     = vr * std::cos(delta) + vi * std::sin(delta) + id * Xqpp_ - iq * Ra_;
      ScalarT psiqpp = -vd / (ONE<RealT> + omega);
      ScalarT psidpp = vq / (ONE<RealT> + omega);
      ScalarT Te     = (psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id;
      ScalarT psiqp  = -(-(Xqp_ - Xl_) * iq + psiqpp * (Xqp_ - Xl_) / (Xqpp_ - Xl_))
                      / (ONE<RealT> + (Xqp_ - Xqpp_) / (Xqpp_ - Xl_));
      ScalarT Edp   = psiqp - (Xqp_ - Xl_) * iq;
      ScalarT psidp = -((Xdp_ - Xl_) * id - psidpp * (Xdp_ - Xl_) / (Xdpp_ - Xl_))
                      / (ONE<RealT> + (Xdp_ - Xdpp_) / (Xdpp_ - Xl_));
      ScalarT Eqp = psidp + (Xdp_ - Xl_) * id;

      /* Now we have the state variables, solve for alg. variables */
      ScalarT ksat;
      ScalarT psipp;

      y_[0] = delta; //= 0.55399038;
      y_[1] = omega; // = 0;
      y_[2] = Eqp;   // = 0.995472581;
      y_[3] = psidp; // = 0.971299567;
      y_[4] = psiqp; // = 0.306880069;
      y_[5] = Edp;   // = 0;

      y_[6] = psiqpp = -psiqp * Xq4_ - Edp * Xq5_;
      y_[7] = psidpp = psidp * Xd4_ + Eqp * Xd5_;
      y_[8] = psipp = std::sqrt(psiqpp * psiqpp + psidpp * psidpp);
      y_[9] = ksat = SB_ * ((psipp - SA_) * (psipp - SA_));
      y_[10] = vd = -psiqpp * (ONE<RealT> + omega);
      y_[11] = vq = psidpp * (ONE<RealT> + omega);
      y_[12] = Te = (psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id;
      y_[13]      = id;
      y_[14]      = iq;
      y_[15]      = ir;
      y_[16]      = ii;
      // y_[17] = efd_set_ = Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat;
      y_[17]      = G_ * (vd * std::sin(delta) + vq * std::cos(delta))
               - B_ * (vd * -std::cos(delta) + vq * std::sin(delta)); /* inort, real */
      y_[18] = B_ * (vd * std::sin(delta) + vq * std::cos(delta))
               + G_ * (vd * -std::cos(delta) + vq * std::sin(delta)); /* inort, imag */

      // Set Setpoint mechanical power, which may or may not be used
      pmech_set_ = Te;
      if (signals_.template isAttached<GenrouExternalVariables::PM>())
      {
        signals_.template writeExternalVariable<GenrouExternalVariables::PM>(Te);
      }

      // Set Efield signal (may or may not exist)
      efd_set_ = Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat;
      if (signals_.template isAttached<GenrouExternalVariables::EFD>())
      {
        signals_.template writeExternalVariable<GenrouExternalVariables::EFD>(efd_set_);
      }

      for (IdxT i = 0; i < size_; ++i)
      {
        yp_[static_cast<size_t>(i)] = 0.0;
      }

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[static_cast<size_t>(i)] = i < 6;
      }
      return 0;
    }

    /**
     * @brief Set absolute tolerance
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::setAbsoluteTolerance()
    {
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline int Genrou<ScalarT, IdxT>::evaluateInternalResidual(
        ScalarT* y,
        ScalarT* yp,
        ScalarT* wb,
        ScalarT* ws,
        ScalarT* f)
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
      ScalarT inr    = y[17];
      ScalarT ini    = y[18];

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
      ScalarT pmech = ws[0];
      ScalarT efd   = ws[1];

      /* 6 Genrou differential equations */
      f[0] = delta_dot - omega * (TWO<RealT> * M_PI * 60.0);
      f[1] = omega_dot - (ONE<RealT> / (TWO<RealT> * H_)) * ((pmech - D_ * omega) / (ONE<RealT> + omega) - telec);
      f[2] = Eqp_dot - (ONE<RealT> / Tdop_) * (efd - (Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat));
      f[3] = psidp_dot - (ONE<RealT> / Tdopp_) * (Eqp - psidp - Xd2_ * id);
      f[4] = psiqp_dot - (ONE<RealT> / Tqopp_) * (Edp - psiqp + Xq2_ * iq);
      f[5] = Edp_dot - (ONE<RealT> / Tqop_) * (-Edp + Xqd_ * psiqpp * ksat + Xq1_ * (iq - Xq3_ * (Edp + iq * Xq2_ - psiqp)));

      /* 11 Genrou algebraic equations */
      f[6]  = psiqpp - (-psiqp * Xq4_ - Edp * Xq5_);
      f[7]  = psidpp - (psidp * Xd4_ + Eqp * Xd5_);
      f[8]  = psipp - std::sqrt((psidpp * psidpp) + (psiqpp * psiqpp));
      f[9]  = ksat - SB_ * ((psipp - SA_) * (psipp - SA_));
      f[10] = vd + psiqpp * (ONE<RealT> + omega);
      f[11] = vq - psidpp * (ONE<RealT> + omega);
      f[12] = telec - ((psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id);
      f[13] = id - (ir * std::sin(delta) - ii * std::cos(delta));
      f[14] = iq - (ir * std::cos(delta) + ii * std::sin(delta));
      f[15] = ir + G_ * vr - B_ * vi - inr;
      f[16] = ii + B_ * vr + G_ * vi - ini;

      /* 2 Genrou current source definitions */
      f[17] = inr - (G_ * (std::sin(delta) * vd + std::cos(delta) * vq) - B_ * (-std::cos(delta) * vd + std::sin(delta) * vq));
      f[18] = ini - (B_ * (std::sin(delta) * vd + std::cos(delta) * vq) + G_ * (-std::cos(delta) * vd + std::sin(delta) * vq));

      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline int Genrou<ScalarT, IdxT>::evaluateBusResidual(
        ScalarT*                  y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      ScalarT inr = y[17];
      ScalarT ini = y[18];
      ScalarT vr  = wb[0];
      ScalarT vi  = wb[1];

      h[0] = (inr - vr * G_ + vi * B_) * mva_base_ / mva_system_base_;
      h[1] = (ini - vr * B_ - vi * G_) * mva_base_ / mva_system_base_;

      return 0;
    }

    /**
     * \brief Residual evaluation and contribution to the connected bus
     *
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::evaluateResidual()
    {
      // Mechanical Power
      ws_[0] = pmech_set_;
      if (signals_.template isAttached<GenrouExternalVariables::PM>())
      {
        ws_[0]         = signals_.template readExternalVariable<GenrouExternalVariables::PM>();
        ws_indices_[0] = signals_.template readExternalVariableIndex<GenrouExternalVariables::PM>();
      }

      // Exciter Efield
      ws_[1] = efd_set_;
      if (signals_.template isAttached<GenrouExternalVariables::EFD>())
      {
        ws_[1]         = signals_.template readExternalVariable<GenrouExternalVariables::EFD>();
        ws_indices_[1] = signals_.template readExternalVariableIndex<GenrouExternalVariables::EFD>();
      }

      // Bus voltages
      wb_[0] = Vr();
      wb_[1] = Vi();

      // Residual evaluation
      evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
      evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());

      // Genrou contribution to bus algebraic equations
      Ir() += h_[0];
      Ii() += h_[1];

      return 0;
    }

    /**
     * @brief Access generator relative speed
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    ScalarT Genrou<ScalarT, IdxT>::getSpeed()
    {
      return y_[1];
    }

    template <class ScalarT, typename IdxT>
    ScalarT Genrou<ScalarT, IdxT>::getTorque()
    {
      return y_[12];
    }

    template <class ScalarT, typename IdxT>
    void Genrou<ScalarT, IdxT>::setDerivedParams()
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
    }
  } // namespace PhasorDynamics
} // namespace GridKit

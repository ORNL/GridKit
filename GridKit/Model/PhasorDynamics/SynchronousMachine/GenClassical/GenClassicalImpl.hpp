/**
 * @file GenClassicalImpl.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a Classical generator model.
 *
 *
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassical.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassicalData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a classical generator model
     */
    template <class ScalarT, typename IdxT>
    GenClassical<ScalarT, IdxT>::GenClassical(bus_type* bus, int unit_id)
      : bus_(bus),
        bus_id_(0),
        unit_id_(unit_id),
        p0_(0.0),
        q0_(0.0),
        H_(3.0),
        D_(0.0),
        Ra_(0.0),
        Xdp_(0.5),
        mva_base_(100.)
    {
      size_ = 5;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a classical generator model
     */
    template <class ScalarT, typename IdxT>
    GenClassical<ScalarT, IdxT>::GenClassical(bus_type* bus,
                                              int       unit_id,
                                              RealT     p0,
                                              RealT     q0,
                                              RealT     H,
                                              RealT     D,
                                              RealT     Ra,
                                              RealT     Xdp)
      : bus_(bus),
        bus_id_(0),
        unit_id_(unit_id),
        p0_(p0),
        q0_(q0),
        H_(H),
        D_(D),
        Ra_(Ra),
        Xdp_(Xdp),
        mva_base_(100.)
    {
      size_ = 5;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a classical generator model
     */
    template <class ScalarT, typename IdxT>
    GenClassical<ScalarT, IdxT>::GenClassical(bus_type* bus, const DataT& data)
      : bus_(bus),
        unit_id_(1),
        monitor_(std::make_unique<MonitorT>(data))
    {
      if (data.parameters.contains(DataT::Parameters::p0))
      {
        p0_ = std::get<RealT>(data.parameters.at(DataT::Parameters::p0));
      }

      if (data.parameters.contains(DataT::Parameters::q0))
      {
        q0_ = std::get<RealT>(data.parameters.at(DataT::Parameters::q0));
      }

      if (data.parameters.contains(DataT::Parameters::H))
      {
        H_ = std::get<RealT>(data.parameters.at(DataT::Parameters::H));
      }

      if (data.parameters.contains(DataT::Parameters::D))
      {
        D_ = std::get<RealT>(data.parameters.at(DataT::Parameters::D));
      }

      if (data.parameters.contains(DataT::Parameters::Ra))
      {
        Ra_ = std::get<RealT>(data.parameters.at(DataT::Parameters::Ra));
      }

      if (data.parameters.contains(DataT::Parameters::Xdp))
      {
        Xdp_ = std::get<RealT>(data.parameters.at(DataT::Parameters::Xdp));
      }

      if (data.parameters.contains(DataT::Parameters::mva_base))
      {
        mva_base_ = std::get<RealT>(data.parameters.at(DataT::Parameters::mva_base));
      }

      if (data.ports.contains(DataT::Ports::bus))
      {
        bus_id_ = data.ports.at(DataT::Ports::bus);
      }

      initializeMonitor();

      size_ = 5;
      setDerivedParams();
    }

    template <class ScalarT, typename IdxT>
    GenClassical<ScalarT, IdxT>::~GenClassical()
    {
    }

    template <class ScalarT, typename IdxT>
    const Model::VariableMonitorBase* GenClassical<ScalarT, IdxT>::getMonitor() const
    {
      return monitor_.get();
    }

    template <class ScalarT, typename IdxT>
    void GenClassical<ScalarT, IdxT>::initializeMonitor()
    {
      using Variable = typename DataT::MonitorableVariables;
      monitor_->set(Variable::ir, [this]
                    { return y_[3]; });
      monitor_->set(Variable::ii, [this]
                    { return y_[4]; });
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
    int GenClassical<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::allocate()
    {
      // Resize component model data
      auto size = static_cast<size_t>(size_);
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      absTol_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Resize coupling data
      wb_.resize(2);
      h_.resize(2);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      return 0;
    }

    /**
     * Initialization of the generator model
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::initialize()
    {
      ScalarT vr    = Vr();
      ScalarT vi    = Vi();
      ScalarT p     = static_cast<ScalarT>(p0_) * mva_system_base_ / mva_base_;
      ScalarT q     = static_cast<ScalarT>(q0_) * mva_system_base_ / mva_base_;
      ScalarT vm2   = vr * vr + vi * vi;
      ScalarT ir    = (p * vr + q * vi) / vm2;
      ScalarT ii    = (p * vi - q * vr) / vm2;
      ScalarT Er    = Ra_ * ir - Xdp_ * ii + vr;
      ScalarT Ei    = Ra_ * ii + Xdp_ * ir + vi;
      ScalarT delta = std::atan2(Ei, Er);
      ScalarT omega = static_cast<ScalarT>(0.0);
      ScalarT Ep    = std::sqrt(Er * Er + Ei * Ei);
      ScalarT Te    = G_ * Ep * Ep - Ep * ((G_ * vr - B_ * vi) * std::cos(delta) + (B_ * vr + G_ * vi) * std::sin(delta));

      y_[0]      = delta;
      y_[1]      = omega;
      y_[2]      = Te;
      y_[3]      = ir;
      y_[4]      = ii;
      pmech_set_ = Te;
      ep_set_    = Ep;

      for (size_t i = 0; i < static_cast<size_t>(size_); ++i)
        yp_[i] = 0.0;

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[static_cast<size_t>(i)] = i < 2;
      }
      return 0;
    }

    /**
     * @brief Set absolute tolerance
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::setAbsoluteTolerance()
    {
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int GenClassical<ScalarT, IdxT>::evaluateInternalResidual(
        ScalarT* y,
        ScalarT* yp,
        ScalarT* wb,
        ScalarT* f)
    {
      // Set variable aliases for better readability.
      const ScalarT delta = y[0];
      const ScalarT omega = y[1];
      const ScalarT telec = y[2];
      const ScalarT ir    = y[3];
      const ScalarT ii    = y[4];
      const ScalarT pmech = pmech_set_; /* Later optionally acquire from governor */
      const ScalarT ep    = ep_set_;    /* Later optionally acquire from exciter */

      // Set derivative aliases for better readability
      const ScalarT delta_dot = yp[0];
      const ScalarT omega_dot = yp[1];

      // Set coupling variable aliases
      const ScalarT vr = wb[0];
      const ScalarT vi = wb[1];

      // GenClassical differential equations
      f[0] = delta_dot - omega * (TWO<RealT> * M_PI * 60.0);
      f[1] = omega_dot - (ONE<RealT> / (TWO<RealT> * H_)) * ((pmech - D_ * omega) / (ONE<RealT> + omega) - telec);

      // GenClassical algebraic equations
      f[2] = telec - (G_ * ep * ep - ep * ((G_ * vr - B_ * vi) * std::cos(delta) + (B_ * vr + G_ * vi) * std::sin(delta)));

      f[3] = ir + G_ * vr - B_ * vi - ep * (G_ * std::cos(delta) - B_ * std::sin(delta));
      f[4] = ii + B_ * vr + G_ * vi - ep * (B_ * std::cos(delta) + G_ * std::sin(delta));

      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int GenClassical<ScalarT, IdxT>::evaluateBusResidual(
        ScalarT*                  y,
        [[maybe_unused]] ScalarT* yp,
        [[maybe_unused]] ScalarT* wb,
        ScalarT*                  h)
    {
      const ScalarT ir = y[3];
      const ScalarT ii = y[4];
      h[0]             = ir * mva_base_ / mva_system_base_;
      h[1]             = ii * mva_base_ / mva_system_base_;

      return 0;
    }

    /**
     * \brief Residual for the generator model.
     *
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::evaluateResidual()
    {
      wb_[0] = Vr();
      wb_[1] = Vi();

      evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), f_.data());
      evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());

      Ir() += h_[0];
      Ii() += h_[1];

      return 0;
    }

    template <class ScalarT, typename IdxT>
    void GenClassical<ScalarT, IdxT>::setDerivedParams()
    {
      G_ = Ra_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
      B_ = -Xdp_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
    }

  } // namespace PhasorDynamics
} // namespace GridKit

/**
 * @file GenClassicalImpl.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a Classical generator model.
 *
 *
 */

#pragma once

#include <algorithm>
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
    template <typename scalar_type, typename index_type>
    GenClassical<scalar_type, index_type>::GenClassical(BusT* bus)
      : bus_(bus),
        bus_id_(0),
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
    template <typename scalar_type, typename index_type>
    GenClassical<scalar_type, index_type>::GenClassical(BusT* bus,
                                                        RealT p0,
                                                        RealT q0,
                                                        RealT H,
                                                        RealT D,
                                                        RealT Ra,
                                                        RealT Xdp)
      : bus_(bus),
        bus_id_(0),
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
    template <typename scalar_type, typename index_type>
    GenClassical<scalar_type, index_type>::GenClassical(BusT* bus, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
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

      if (data.parameters.contains(Parameter::Xdp))
      {
        Xdp_ = std::get<RealT>(data.parameters.at(Parameter::Xdp));
      }

      if (data.parameters.contains(Parameter::mva))
      {
        mva_base_ = std::get<RealT>(data.parameters.at(Parameter::mva));
      }

      if (data.buses.contains(Buses::bus))
      {
        bus_id_ = data.buses.at(Buses::bus);
      }

      initializeMonitor();

      size_ = 5;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    GenClassical<scalar_type, index_type>::~GenClassical()
    {
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* GenClassical<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    GenClassical<scalar_type, index_type>::ScalarT GenClassical<scalar_type, index_type>::toMachineBase(ScalarT value) const
    {
      return value * va_system_base_ / va_machine_base_;
    }

    template <typename scalar_type, typename index_type>
    GenClassical<scalar_type, index_type>::ScalarT GenClassical<scalar_type, index_type>::toSystemBase(ScalarT value) const
    {
      return value / toMachineBase(static_cast<ScalarT>(ONE<RealT>));
    }

    template <typename scalar_type, typename index_type>
    void GenClassical<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      monitor_->set(Variable::ir, [this]
                    { return toSystemBase(y_.getData()[3]); });
      monitor_->set(Variable::ii, [this]
                    { return toSystemBase(y_.getData()[4]); });
      monitor_->set(Variable::p, [this]
                    { return toSystemBase(Vr() * y_.getData()[3] + Vi() * y_.getData()[4]); });
      monitor_->set(Variable::q, [this]
                    { return toSystemBase(Vi() * y_.getData()[3] - Vr() * y_.getData()[4]); });
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
    int GenClassical<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::allocate()
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

      // Resize coupling data
      this->allocateExternalVectors(static_cast<IdxT>(GenClassicalExternalVariables::MAXIMUM));
      f_ext_.resize(2);
      residual_indices_ext_.assign(2, INVALID_INDEX<IdxT>);

      allocated_ = true;
      return 0;
    }

    /**
     * Initialization of the generator model
     */
    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::initialize()
    {
      ScalarT vr    = Vr();
      ScalarT vi    = Vi();
      ScalarT p     = toMachineBase(static_cast<ScalarT>(p0_));
      ScalarT q     = toMachineBase(static_cast<ScalarT>(q0_));
      ScalarT vm2   = vr * vr + vi * vi;
      ScalarT ir    = (p * vr + q * vi) / vm2;
      ScalarT ii    = (p * vi - q * vr) / vm2;
      ScalarT Er    = Ra_ * ir - Xdp_ * ii + vr;
      ScalarT Ei    = Ra_ * ii + Xdp_ * ir + vi;
      ScalarT delta = std::atan2(Ei, Er);
      ScalarT omega = static_cast<ScalarT>(0.0);
      ScalarT Ep    = std::sqrt(Er * Er + Ei * Ei);
      ScalarT Te    = G_ * Ep * Ep - Ep * ((G_ * vr - B_ * vi) * std::cos(delta) + (B_ * vr + G_ * vi) * std::sin(delta));

      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]       = delta;
      y[1]       = omega;
      y[2]       = Te;
      y[3]       = ir;
      y[4]       = ii;
      pmech_set_ = Te;
      ep_set_    = Ep;

      for (size_t i = 0; i < static_cast<size_t>(size_); ++i)
        yp[i] = 0.0;

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[static_cast<size_t>(i)] = i < 2;
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
    int GenClassical<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int GenClassical<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y,
        const ScalarT* yp,
        const ScalarT* y_ext,
        ScalarT*       f)
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
      const ScalarT vr = y_ext[0];
      const ScalarT vi = y_ext[1];

      static constexpr auto pi = std::numbers::pi_v<RealT>;

      // GenClassical differential equations
      f[0] = delta_dot - omega * (TWO<RealT> * pi * freq_system_base_);
      f[1] = omega_dot - H_inv_ * ((pmech - D_ * omega) / (ONE<RealT> + omega) - telec);

      // GenClassical algebraic equations
      f[2] = telec - (G_ * ep * ep - ep * ((G_ * vr - B_ * vi) * std::cos(delta) + (B_ * vr + G_ * vi) * std::sin(delta)));

      f[3] = ir + G_ * vr - B_ * vi - ep * (G_ * std::cos(delta) - B_ * std::sin(delta));
      f[4] = ii + B_ * vr + G_ * vi - ep * (B_ * std::cos(delta) + G_ * std::sin(delta));

      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int GenClassical<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT ir = y[3];
      const ScalarT ii = y[4];
      f_ext[0]         = toSystemBase(ir);
      f_ext[1]         = toSystemBase(ii);

      return 0;
    }

    /**
     * @brief Gather external variables and index maps.
     *
     */
    template <typename scalar_type, typename index_type>
    void GenClassical<scalar_type, index_type>::gatherExternalVariables()
    {
      y_ext_[0] = Vr();
      y_ext_[1] = Vi();
      if (bus_->size() > 0)
      {
        variable_indices_ext_[0] = bus_->getVariableIndex(0);
        variable_indices_ext_[1] = bus_->getVariableIndex(1);
        residual_indices_ext_[0] = bus_->getResidualIndex(0);
        residual_indices_ext_[1] = bus_->getResidualIndex(1);
      }
    }

    /**
     * \brief Internal residual for the generator model.
     *
     */
    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::evaluateInternalResidual()
    {
      gatherExternalVariables();

      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      evaluateInternalResidual(y, yp, y_ext_.data(), f);
      f_.setDataUpdated();

      return 0;
    }

    /**
     * \brief External residual contributions to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::evaluateExternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), f_ext_.data());

      Ir() += f_ext_[0];
      Ii() += f_ext_[1];

      if (bus_->size() > 0)
      {
        bus_->getResidual().setDataUpdated();
      }

      return 0;
    }

    /**
     * \brief Residual for the generator model.
     *
     */
    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
    }

    template <typename scalar_type, typename index_type>
    void GenClassical<scalar_type, index_type>::setDerivedParams()
    {
      if (H_ == ZERO<RealT> && D_ == ZERO<RealT>)
      {
        H_inv_ = ZERO<RealT>;
      }
      else
      {
        H_     = std::max(H_, static_cast<RealT>(0.1));
        H_inv_ = ONE<RealT> / (TWO<RealT> * H_);
      }

      Xdp_ = std::clamp(Xdp_, static_cast<RealT>(1.0e-5), static_cast<RealT>(999.0));

      G_               = Ra_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
      B_               = -Xdp_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
      va_machine_base_ = mva_base_ * static_cast<RealT>(1.0e6);
    }

  } // namespace PhasorDynamics
} // namespace GridKit

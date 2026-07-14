#pragma once

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Constructor for a bus fault
     *
     * Model sizes:
     * - Number of equations = 0
     * - Number of independent variables = 0
     */
    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>::BusFault(BusT* bus)
      : bus_(bus), R_(0), X_(0.01), status_(0), bus_id_(0)
    {
      (void) bus_id_;
      size_ = 2;
      setDerivedParams();
    }

    /**
     * @brief Construct a new BusFault
     *
     * @param bus1 - pointer to bus-1
     * @param bus2 - pointer to bus-2
     * @param R - line series resistance
     * @param X - line series reactance
     * @param G - line shunt conductance
     * @param B - line shunt charging
     */
    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>::BusFault(BusT* bus, RealT R, RealT X, int status)
      : bus_(bus), R_(R), X_(X), status_(status), bus_id_(0)
    {
      size_ = 2;
      setDerivedParams();
    }

    /**
     * @brief Construct a new BusFault
     *
     * @param bus1 - pointer to bus-1
     * @param bus2 - pointer to bus-2
     */
    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>::BusFault(BusT* bus, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      using Parameter = typename ModelDataT::Parameters;
      using Buses     = typename ModelDataT::Buses;

      if (data.parameters.contains(Parameter::R))
      {
        R_ = std::get<RealT>(data.parameters.at(Parameter::R));
      }

      if (data.parameters.contains(Parameter::X))
      {
        X_ = std::get<RealT>(data.parameters.at(Parameter::X));
      }

      if (data.parameters.contains(Parameter::state0))
      {
        status_ = std::get<bool>(data.parameters.at(Parameter::state0));
      }

      if (data.buses.contains(Buses::bus))
      {
        bus_id_ = data.buses.at(Buses::bus);
      }

      using Variable = typename ModelDataT::MonitorableVariables;
      monitor_->set(Variable::state, [this]
                    { return status_; });
      monitor_->set(Variable::ir, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::ii, [this]
                    { return y_.getData()[1]; });

      size_ = 2;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>::~BusFault()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::allocate()
    {
      if (!this->allocated_)
      {
        this->allocateVectors(this->size_);
      }
      // std::cout << "Allocate BusFault..." << std::endl;
      auto size = static_cast<std::size_t>(size_);

      tag_.resize(size);

      this->variable_indices_.resize(size);
      this->residual_indices_.resize(size);

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
     * Initialization of the branch model
     *
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::initialize()
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      if (status_)
      {
        ScalarT vr = Vr();
        ScalarT vi = Vi();
        ScalarT ir = -(vr * G_ - vi * B_);
        ScalarT ii = -(vr * B_ + vi * G_);
        y[0]       = ir;
        y[1]       = ii;
      }
      else
      {
        y[0] = 0.0;
        y[1] = 0.0;
      }

      yp[0] = 0.0;
      yp[1] = 0.0;

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = false;
      tag_[1] = false;

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
    template <class scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int BusFault<scalar_type, index_type>::evaluateBusResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* wb,
        ScalarT*                        h)
    {
      const ScalarT Ir = y[0];
      const ScalarT Ii = y[1];
      h[0]             = Ir;
      h[1]             = Ii;

      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int BusFault<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        f)
    {
      const ScalarT Vr = wb[0];
      const ScalarT Vi = wb[1];
      const ScalarT Ir = y[0];
      const ScalarT Ii = y[1];
      f[0]             = Ir + Vr * G_ - Vi * B_;
      f[1]             = Ii + Vr * B_ + Vi * G_;

      return 0;
    }

    /**
     * \brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     *
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::evaluateResidual()
    {
      if (status_)
      {
        wb_[0]         = Vr();
        wb_[1]         = Vi();
        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb_.data(), f);
        evaluateBusResidual(y, yp, wb_.data(), h_.data());
        Ir() += h_[0];
        Ii() += h_[1];
      }
      else
      {
        wb_[0]         = 0.0;
        wb_[1]         = 0.0;
        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb_.data(), f);
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* BusFault<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename scalar_type, typename index_type>
    void BusFault<scalar_type, index_type>::setDerivedParams()
    {
      B_ = -X_ / (X_ * X_ + R_ * R_);
      G_ = R_ / (X_ * X_ + R_ * R_);
    }
  } // namespace PhasorDynamics
} // namespace GridKit

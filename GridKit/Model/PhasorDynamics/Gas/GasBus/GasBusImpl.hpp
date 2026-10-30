/**
 * @file GasBusImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the gas bus model.
 */
#pragma once

#include <cmath>

#include <GridKit/Model/PhasorDynamics/Gas/GasBus/GasBus.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Construct a gas bus
       *
       * Model size:
       * - Number of equations = 1
       * - Number of internal variables = 1
       *
       * @param data - gas bus data with the initial pressure and monitor selections
       */
      template <typename scalar_type, typename index_type>
      GasBus<scalar_type, index_type>::GasBus(const ModelDataT& data)
        : p0_(data.p0),
          monitor_(std::make_unique<MonitorT>("GasBus", data.name, data.monitored_variables))
      {
        bus_id_ = data.bus_id;
        size_   = 1;

        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::p, [this]
                      { return p(); });
      }

      template <typename scalar_type, typename index_type>
      GasBus<scalar_type, index_type>::~GasBus()
      {
      }

      template <typename scalar_type, typename index_type>
      int GasBus<scalar_type, index_type>::allocate()
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

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks GasBus<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("GasBus");

        checks.check(std::isfinite(p0_) && p0_ > ZERO<RealT>,
                     "initial pressure must be finite and positive");

        return checks;
      }

      /**
       * @brief Set the bus pressure to its initial value
       */
      template <typename scalar_type, typename index_type>
      int GasBus<scalar_type, index_type>::initialize()
      {
        if (!this->verifyAndLog())
        {
          return 1;
        }

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        y[0]  = p0_;
        yp[0] = ZERO<RealT>;

        y_.setDataUpdated();
        yp_.setDataUpdated();

        if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
        {
          this->initializeDependencyTrackingVariableNumbers();
        }

        return 0;
      }

      /**
       * @brief The pressure is differential: connected gas branches add its
       * derivative to the mass-flow balance.
       */
      template <typename scalar_type, typename index_type>
      int GasBus<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[0] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GasBus<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Reset the mass-flow balance
       *
       * @warning This implementation assumes the gas bus residual is evaluated
       * _before_ the residuals of the devices connected to it.
       */
      template <typename scalar_type, typename index_type>
      int GasBus<scalar_type, index_type>::evaluateResidual()
      {
        q() = ZERO<RealT>;
        f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* GasBus<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit

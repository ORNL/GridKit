/**
 * @file GasBusInfiniteImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the fixed-pressure gas bus model.
 */
#pragma once

#include <cmath>

#include <GridKit/Model/PhasorDynamics/Gas/GasBusInfinite/GasBusInfinite.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Construct a fixed-pressure gas bus
       *
       * Model size:
       * - Number of equations = 0
       * - Number of internal variables = 0
       *
       * @param data - gas bus data with the pressure and monitor selections
       */
      template <typename scalar_type, typename index_type>
      GasBusInfinite<scalar_type, index_type>::GasBusInfinite(const ModelDataT& data)
        : p0_(data.p0),
          p_(data.p0),
          monitor_(std::make_unique<MonitorT>("GasBusInfinite", data.name, data.monitored_variables))
      {
        bus_id_ = data.bus_id;
        size_   = 0;

        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::p, [this]
                      { return p(); });
      }

      template <typename scalar_type, typename index_type>
      GasBusInfinite<scalar_type, index_type>::~GasBusInfinite()
      {
      }

      template <typename scalar_type, typename index_type>
      int GasBusInfinite<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        variable_indices_.resize(size);
        residual_indices_.resize(size);

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks GasBusInfinite<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("GasBusInfinite");

        checks.check(std::isfinite(p0_) && p0_ > ZERO<RealT>,
                     "pressure must be finite and positive");

        return checks;
      }

      /**
       * @brief Restore the case-file pressure
       */
      template <typename scalar_type, typename index_type>
      int GasBusInfinite<scalar_type, index_type>::initialize()
      {
        if (!this->verifyAndLog())
        {
          return 1;
        }

        p_ = p0_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GasBusInfinite<scalar_type, index_type>::tagDifferentiable()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GasBusInfinite<scalar_type, index_type>::setAbsoluteTolerance(RealT)
      {
        return 0;
      }

      /**
       * @brief Discard the mass flows of the previous evaluation
       */
      template <typename scalar_type, typename index_type>
      int GasBusInfinite<scalar_type, index_type>::evaluateResidual()
      {
        q_ = ZERO<RealT>;
        return 0;
      }

      /**
       * @brief The bus has no variables and no Jacobian entries
       */
      template <typename scalar_type, typename index_type>
      int GasBusInfinite<scalar_type, index_type>::evaluateJacobian()
      {
        return this->constructCoo();
      }

      /**
       * @brief Change the fixed pressure, for example at an event
       *
       * @param[in] pressure - new bus pressure [Pa]
       */
      template <typename scalar_type, typename index_type>
      int GasBusInfinite<scalar_type, index_type>::setPressure(RealT pressure)
      {
        if (!std::isfinite(pressure) || pressure <= ZERO<RealT>)
        {
          Log::error() << "GasBusInfinite: pressure must be finite and positive\n";
          return 1;
        }

        p_ = pressure;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* GasBusInfinite<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit

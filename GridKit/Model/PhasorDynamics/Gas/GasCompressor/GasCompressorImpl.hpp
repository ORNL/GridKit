/**
 * @file GasCompressorImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the gas compressor model.
 */
#pragma once

#include <cmath>
#include <limits>

#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasCompressor/GasCompressor.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Construct a gas compressor
       *
       * Model size:
       * - Number of equations = 1
       * - Number of internal variables = 1
       *
       * @param bus1 - suction bus
       * @param bus2 - discharge bus
       * @param data - gas compressor parameters and monitored-variable selections
       */
      template <typename scalar_type, typename index_type>
      GasCompressor<scalar_type, index_type>::GasCompressor(GasBusT* bus1, GasBusT* bus2, const ModelDataT& data)
        : bus1_(bus1),
          bus2_(bus2),
          monitor_(std::make_unique<MonitorT>(data))
      {
        size_ = 1;
        initializeParameters(data);

        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::q, [this]
                      { return y_.getData()[0]; });
      }

      template <typename scalar_type, typename index_type>
      GasCompressor<scalar_type, index_type>::~GasCompressor()
      {
      }

      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void GasCompressor<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;

        const Model::ParameterReader reader(data, "GasCompressor");

        ratio_ = reader.readReal(Parameter::ratio).value_or(ratio_);
      }

      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.resize(size);

        variable_indices_.resize(size);
        residual_indices_.resize(size);
        this->setVariableIndex(0, 0);
        this->setResidualIndex(0, 0);

        wbp_.resize(2);
        h_.resize(2);

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks GasCompressor<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("GasCompressor");

        checks.check(bus1_ != nullptr, "bus1 pointer is null");
        checks.check(bus2_ != nullptr, "bus2 pointer is null");
        if (bus1_ != nullptr && bus2_ != nullptr)
        {
          checks.check(bus1_->size() > 0 || bus2_->size() > 0,
                       "both buses have a fixed pressure");
        }

        checks.check(std::isfinite(ratio_) && ratio_ > ZERO<RealT>,
                     "ratio must be finite and positive");

        return checks;
      }

      /**
       * @brief Check the terminal pressures against the ratio; the flow starts
       * at zero for the consistent initial condition solve
       */
      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::initialize()
      {
        if (!this->verifyAndLog())
        {
          return 1;
        }

        // The ratio is enforced through its derivative, which keeps any
        // initial mismatch. Full-precision case pressures agree to a few ulps.
        const RealT   tolerance = 100 * std::numeric_limits<RealT>::epsilon();
        const ScalarT p1        = bus1_->p();
        const ScalarT p2        = bus2_->p();
        if (p1 <= ZERO<RealT> || p2 <= ZERO<RealT> || std::abs(ratio_ * p1 - p2) > tolerance * p2)
        {
          Log::error() << "GasCompressor: terminal-bus pressures must be positive with p2 = ratio * p1\n";
          return 1;
        }

        y_.getData()[0]  = ZERO<RealT>;
        yp_.getData()[0] = ZERO<RealT>;

        y_.setDataUpdated();
        yp_.setDataUpdated();

        if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
        {
          this->initializeDependencyTrackingVariableNumbers();
        }

        return 0;
      }

      /// The mass flow is algebraic
      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[0] = false;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /// Read the bus pressure derivatives; a fixed-pressure bus has none
      template <typename scalar_type, typename index_type>
      void GasCompressor<scalar_type, index_type>::gatherBusValues()
      {
        auto* wbp = wbp_.getData();

        wbp[0] = ZERO<RealT>;
        if (bus1_->size() > 0)
        {
          wbp[0] = bus1_->yp().getData()[0];
        }

        wbp[1] = ZERO<RealT>;
        if (bus2_->size() > 0)
        {
          wbp[1] = bus2_->yp().getData()[0];
        }
      }

      /**
       * @brief Time derivative of the pressure ratio
       */
      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::evaluateInternalResidual(const ScalarT* wbp, ScalarT* f)
      {
        f[0] = ratio_ * wbp[0] - wbp[1];
        return 0;
      }

      /**
       * @brief Mass flow leaving bus 1 and entering bus 2
       */
      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::evaluateBusResidual(const ScalarT* y, ScalarT* h)
      {
        h[0] = -y[0];
        h[1] = y[0];
        return 0;
      }

      /**
       * @brief Evaluate the ratio residual and add the mass flow to the buses
       */
      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::evaluateResidual()
      {
        gatherBusValues();

        auto* h = h_.getData();
        evaluateInternalResidual(wbp_.getData(), f_.getData());
        evaluateBusResidual(y_.getData(), h);

        bus1_->q() += h[0];
        bus2_->q() += h[1];

        if (bus1_->size() > 0)
        {
          bus1_->getResidual().setDataUpdated();
        }
        if (bus2_->size() > 0)
        {
          bus2_->getResidual().setDataUpdated();
        }
        f_.setDataUpdated();

        return 0;
      }

      /// Global row of a bus balance; a fixed-pressure bus has none
      template <typename scalar_type, typename index_type>
      index_type GasCompressor<scalar_type, index_type>::busRow(GasBusT* bus)
      {
        if (bus->size() == 0)
        {
          return INVALID_INDEX<IdxT>;
        }
        return bus->getResidualIndex(0);
      }

      /// Global column of a bus pressure; a fixed-pressure bus has none
      template <typename scalar_type, typename index_type>
      index_type GasCompressor<scalar_type, index_type>::busColumn(GasBusT* bus)
      {
        if (bus->size() == 0)
        {
          return INVALID_INDEX<IdxT>;
        }
        return bus->getVariableIndex(0);
      }

      /// Store one Jacobian entry; entries of a fixed-pressure bus are dropped
      template <typename scalar_type, typename index_type>
      void GasCompressor<scalar_type, index_type>::addEntry(IdxT row, IdxT column, RealT value)
      {
        if (row == INVALID_INDEX<IdxT> || column == INVALID_INDEX<IdxT>)
        {
          return;
        }

        J_rows_buffer_[nnz_] = row;
        J_cols_buffer_[nnz_] = column;
        J_vals_buffer_[nnz_] = value;
        ++nnz_;
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* GasCompressor<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit

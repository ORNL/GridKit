/**
 * @file GasBranchImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the gas branch (pipe) model.
 */
#pragma once

#include <cmath>
#include <numbers>

#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBranch/GasBranch.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Construct a gas branch
       *
       * Model size:
       * - Number of equations = 2N - 1
       * - Number of internal variables = 2N - 1
       *
       * @param bus1 - bus-1 terminal
       * @param bus2 - bus-2 terminal
       * @param data - gas branch parameters and monitored-variable selections
       */
      template <typename scalar_type, typename index_type>
      GasBranch<scalar_type, index_type>::GasBranch(GasBusT* bus1, GasBusT* bus2, const ModelDataT& data)
        : bus1_(bus1),
          bus2_(bus2),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        setDerivedParameters();

        size_ = 2 * N_ - 1;
      }

      template <typename scalar_type, typename index_type>
      GasBranch<scalar_type, index_type>::~GasBranch()
      {
      }

      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void GasBranch<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;

        const Model::ParameterReader reader(data, "GasBranch");

        L_      = reader.readReal(Parameter::L).value_or(L_);
        D_      = reader.readReal(Parameter::D).value_or(D_);
        lambda_ = reader.readReal(Parameter::lambda).value_or(lambda_);
        a_      = reader.readReal(Parameter::a).value_or(a_);
        dxmax_  = reader.readReal(Parameter::dxmax).value_or(dxmax_);
      }

      /**
       * @brief Segment count, segment length, and pipe area
       */
      template <typename scalar_type, typename index_type>
      void GasBranch<scalar_type, index_type>::setDerivedParameters()
      {
        // verify() rejects the parameters that leave N at its minimum
        N_ = 1;
        if (L_ > ZERO<RealT> && dxmax_ > ZERO<RealT>)
        {
          N_ = static_cast<IdxT>(std::ceil(L_ / dxmax_));
        }
        dx_ = L_ / static_cast<RealT>(N_);
        A_  = std::numbers::pi_v<RealT> * D_ * D_ / FOUR<RealT>;
      }

      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::allocate()
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

        wb_.resize(2);
        wbp_.resize(2);
        h_.resize(2);

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks GasBranch<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("GasBranch");

        checks.check(bus1_ != nullptr, "bus1 pointer is null");
        checks.check(bus2_ != nullptr, "bus2 pointer is null");

        checks.check(std::isfinite(L_) && L_ > ZERO<RealT>, "L must be finite and positive");
        checks.check(std::isfinite(D_) && D_ > ZERO<RealT>, "D must be finite and positive");
        checks.check(std::isfinite(lambda_) && lambda_ > ZERO<RealT>,
                     "lambda must be finite and positive");
        checks.check(std::isfinite(a_) && a_ > ZERO<RealT>, "a must be finite and positive");
        checks.check(std::isfinite(dxmax_) && dxmax_ > ZERO<RealT>,
                     "dxmax must be finite and positive");

        return checks;
      }

      /**
       * @brief Initialize to the steady flow between the terminal pressures
       */
      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::initialize()
      {
        if (!this->verifyAndLog())
        {
          return 1;
        }

        const ScalarT p0 = bus1_->p();
        const ScalarT pN = bus2_->p();
        if (p0 <= ZERO<RealT> || pN <= ZERO<RealT>)
        {
          Log::error() << "GasBranch: terminal-bus pressures must be positive\n";
          return 1;
        }

        const ScalarT dp2 = p0 * p0 - pN * pN;
        ScalarT       phi = std::sqrt(D_ * std::abs(dp2) / (lambda_ * a_ * a_ * L_));
        if (dp2 < ZERO<RealT>)
        {
          phi = -phi;
        }

        auto* y  = y_.getData();
        auto* yp = yp_.getData();
        for (IdxT j = 1; j < N_; ++j)
        {
          y[pIndex(j)] = std::sqrt(p0 * p0 - static_cast<RealT>(j) / static_cast<RealT>(N_) * dp2);
        }
        for (IdxT k = 1; k <= N_; ++k)
        {
          y[phiIndex(k)] = phi;
        }
        for (IdxT i = 0; i < size_; ++i)
        {
          yp[i] = ZERO<RealT>;
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();

        if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
        {
          this->initializeDependencyTrackingVariableNumbers();
        }

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::tagDifferentiable()
      {
        for (IdxT i = 0; i < size_; ++i)
        {
          tag_[static_cast<size_t>(i)] = true;
        }
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /// Local index of the pressure at interior node j
      template <typename scalar_type, typename index_type>
      index_type GasBranch<scalar_type, index_type>::pIndex(IdxT j) const
      {
        return j - 1;
      }

      /// Local index of the mass flux in segment k
      template <typename scalar_type, typename index_type>
      index_type GasBranch<scalar_type, index_type>::phiIndex(IdxT k) const
      {
        return N_ + k - 2;
      }

      /// Pressure at node j; nodes 0 and N are the bus pressures
      template <typename scalar_type, typename index_type>
      scalar_type GasBranch<scalar_type, index_type>::nodePressure(IdxT j, const ScalarT* y, const ScalarT* wb) const
      {
        if (j == 0)
        {
          return wb[0];
        }
        if (j == N_)
        {
          return wb[1];
        }
        return y[pIndex(j)];
      }

      /// Read the bus pressures and their derivatives; a fixed-pressure bus has none
      template <typename scalar_type, typename index_type>
      void GasBranch<scalar_type, index_type>::gatherBusValues()
      {
        auto* wb  = wb_.getData();
        auto* wbp = wbp_.getData();

        wb[0] = bus1_->p();
        wb[1] = bus2_->p();

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
       * @brief Node continuity and segment momentum
       */
      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::evaluateInternalResidual(const ScalarT* y,
                                                                       const ScalarT* yp,
                                                                       const ScalarT* wb,
                                                                       ScalarT*       f)
      {
        const RealT a2 = a_ * a_;

        for (IdxT j = 1; j < N_; ++j)
        {
          f[pIndex(j)] = -yp[pIndex(j)] + a2 / dx_ * (y[phiIndex(j)] - y[phiIndex(j + 1)]);
        }

        for (IdxT k = 1; k <= N_; ++k)
        {
          const ScalarT p_from  = nodePressure(k - 1, y, wb);
          const ScalarT p_to    = nodePressure(k, y, wb);
          const ScalarT phi     = y[phiIndex(k)];
          const ScalarT rho_bar = HALF<RealT> * (p_from + p_to) / a2;

          f[phiIndex(k)] = -yp[phiIndex(k)] + (p_from - p_to) / dx_
                           - lambda_ / (TWO<RealT> * D_) * phi * std::abs(phi) / rho_bar;
        }

        return 0;
      }

      /**
       * @brief Terminal mass flows, including the half-segment storage at each end
       */
      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::evaluateBusResidual(const ScalarT* y,
                                                                  const ScalarT* wbp,
                                                                  ScalarT*       h)
      {
        const RealT storage = HALF<RealT> * A_ * dx_ / (a_ * a_);

        h[0] = -A_ * y[phiIndex(1)] - storage * wbp[0];
        h[1] = A_ * y[phiIndex(N_)] - storage * wbp[1];

        return 0;
      }

      /**
       * @brief Evaluate the pipe residual and add the terminal mass flows to the buses
       */
      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::evaluateResidual()
      {
        gatherBusValues();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        auto*       h  = h_.getData();
        evaluateInternalResidual(y, yp, wb_.getData(), f);
        evaluateBusResidual(y, wbp_.getData(), h);

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
      index_type GasBranch<scalar_type, index_type>::busRow(GasBusT* bus)
      {
        if (bus->size() == 0)
        {
          return INVALID_INDEX<IdxT>;
        }
        return bus->getResidualIndex(0);
      }

      /// Global column of a bus pressure; a fixed-pressure bus has none
      template <typename scalar_type, typename index_type>
      index_type GasBranch<scalar_type, index_type>::busColumn(GasBusT* bus)
      {
        if (bus->size() == 0)
        {
          return INVALID_INDEX<IdxT>;
        }
        return bus->getVariableIndex(0);
      }

      /// Global column of the pressure at node j
      template <typename scalar_type, typename index_type>
      index_type GasBranch<scalar_type, index_type>::nodeColumn(IdxT j)
      {
        if (j == 0)
        {
          return busColumn(bus1_);
        }
        if (j == N_)
        {
          return busColumn(bus2_);
        }
        return this->getVariableIndex(pIndex(j));
      }

      /// Store one Jacobian entry; entries of a fixed-pressure bus are dropped
      template <typename scalar_type, typename index_type>
      void GasBranch<scalar_type, index_type>::addEntry(IdxT row, IdxT column, RealT value)
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
      const Model::VariableMonitorBase* GasBranch<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void GasBranch<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::q1, [this]
                      {
                        ScalarT h[2];
                        gatherBusValues();
                        evaluateBusResidual(y_.getData(), wbp_.getData(), h);
                        return h[0]; });
        monitor_->set(Variable::q2, [this]
                      {
                        ScalarT h[2];
                        gatherBusValues();
                        evaluateBusResidual(y_.getData(), wbp_.getData(), h);
                        return h[1]; });
        monitor_->set(Variable::linepack, [this]
                      {
                        gatherBusValues();
                        const auto* y            = y_.getData();
                        const auto* wb           = wb_.getData();
                        ScalarT     pressure_sum = HALF<RealT> * (wb[0] + wb[1]);
                        for (IdxT j = 1; j < N_; ++j)
                        {
                          pressure_sum += y[pIndex(j)];
                        }
                        return A_ * dx_ * pressure_sum / (a_ * a_); });
      }
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit

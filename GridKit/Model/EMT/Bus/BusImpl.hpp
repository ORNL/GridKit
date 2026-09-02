#pragma once

#include <algorithm>
#include <iostream>

#include <GridKit/Model/EMT/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Bus/BusData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase EMT bus
     *
     * System sizes:
     * - Number of equations = 3
     * - Number of independent variables = 3
     */
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus()
    {
      size_ = 3;
    }

    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(ScalarT va0, ScalarT vb0, ScalarT vc0)
      : va0_(va0),
        vb0_(vb0),
        vc0_(vc0)
    {
      size_ = 3;
    }

    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : bus_id_(data.bus_id),
        va0_(data.va0),
        vb0_(data.vb0),
        vc0_(data.vc0),
        monitor_(std::make_unique<MonitorT>(data))
    {
      size_ = 3;
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::~Bus()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local storage and binds the voltage port.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize coupling data
      this->allocateExternalVectors(0, 0);

      // Bind the voltage port to the phase variables and residual rows
      auto* y  = y_.getData();
      auto* yp = yp_.getData();
      auto* f  = f_.getData();
      for (IdxT n = 0; n < size_; ++n)
      {
        v_port_.nodes[static_cast<size_t>(n)].set(&y[n],
                                                  &yp[n],
                                                  &f[n],
                                                  &(this->getVariableIndex(n)),
                                                  &(this->getResidualIndex(n)));
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::verify() const
    {
      return 0;
    }

    /**
     * Initialization of the bus model
     *
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize()
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = va0_;
      y[1]  = vb0_;
      y[2]  = vc0_;
      yp[0] = 0.0;
      yp[1] = 0.0;
      yp[2] = 0.0;

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    /**
     * \brief Identify differential variables.
     *
     * A phase voltage is differential when a connected component contributes
     * a voltage derivative to the current-balance residual row, marked on the
     * port node during the connected component's allocation.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = v_port_.nodes[0].hasDerivativeCoupling();
      tag_[1] = v_port_.nodes[1].hasDerivativeCoupling();
      tag_[2] = v_port_.nodes[2].hasDerivativeCoupling();

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
    int Bus<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Internal residual
     *
     * The current-balance rows start at zero; connected components accumulate
     * their injections during the external residual phase.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Bus<scalar_type, index_type>::evaluateInternalResidual(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      f[0] = 0.0;
      f[1] = 0.0;
      f[2] = 0.0;

      return 0;
    }

    /**
     * @brief External residual
     *
     * The bus owns no external residual rows.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Bus<scalar_type, index_type>::evaluateExternalResidual(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        [[maybe_unused]] ScalarT*       f_ext)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateInternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      evaluateInternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f);
      f_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateExternalResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Bus<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Bus<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::va, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::vb, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::vc, [this]
                    { return y_.getData()[2]; });
    }

  } // namespace EMT
} // namespace GridKit

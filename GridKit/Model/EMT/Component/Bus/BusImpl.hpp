#pragma once

#include <algorithm>
#include <iostream>

#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Bus/BusData.hpp>
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
      for (auto& voltage : v_port_)
        voltage.claimProducer();
    }

    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      size_ = 3;
      for (auto& voltage : v_port_)
        voltage.claimProducer();
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
      this->allocateExternalVectors(3, 0);
      for (IdxT n = 0; n < 3; ++n)
        this->setExternalVariableSignal(n, currents_[static_cast<size_t>(n)]);

      // Bind the voltage port to the phase variables and residual rows
      auto* y  = y_.getData();
      auto* yp = yp_.getData();
      auto* f  = f_.getData();
      for (IdxT n = 0; n < size_; ++n)
      {
        v_port_[static_cast<size_t>(n)].set(&y[n],
                                            &yp[n],
                                            &f[n],
                                            &(this->getVariableIndex(n)),
                                            &(this->getResidualIndex(n)));
        if (auto* output = outputs_[static_cast<size_t>(n)])
          output->set(&y[n], &yp[n], &f[n], &this->getVariableIndex(n), &this->getResidualIndex(n));
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
      for (const auto* current : currents_)
        if (current != nullptr && !current->linked())
          return 1;
      return 0;
    }

    /**
     * Initialization of the bus model
     *
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = this->outputValue(outputs, Outputs::va, ZERO<RealT>);
      y[1]  = this->outputValue(outputs, Outputs::vb, ZERO<RealT>);
      y[2]  = this->outputValue(outputs, Outputs::vc, ZERO<RealT>);
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
     * port signal during the connected component's allocation.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::tagDifferentiable()
    {
      for (size_t phase = 0; phase < 3; ++phase)
        tag_[phase] = v_port_[phase].hasDerivativeCoupling()
                      || (outputs_[phase] && outputs_[phase]->hasDerivativeCoupling());

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
      f[0] = y_ext[0];
      f[1] = y_ext[1];
      f[2] = y_ext[2];

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
      this->gatherExternalVariables();
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
    int Bus<scalar_type, index_type>::evaluateJacobian()
    {
      std::array<typename SignalT::GradientT, 3> gradients;
      size_t                                     entries = 0;
      for (size_t phase = 0; phase < 3; ++phase)
      {
        if (currents_[phase])
          currents_[phase]->appendGradient(gradients[phase]);
        entries += gradients[phase].size();
      }
      this->resetJacobianStructure();
      if (entries > jacobian_capacity_)
      {
        delete[] J_rows_buffer_;
        delete[] J_cols_buffer_;
        delete[] J_vals_buffer_;
        J_rows_buffer_     = new IdxT[entries];
        J_cols_buffer_     = new IdxT[entries];
        J_vals_buffer_     = new RealT[entries];
        jacobian_capacity_ = entries;
      }
      nnz_ = 0;
      for (size_t phase = 0; phase < 3; ++phase)
        for (const auto& [column, value] : gradients[phase])
        {
          J_rows_buffer_[nnz_]   = this->getResidualIndex(static_cast<IdxT>(phase));
          J_cols_buffer_[nnz_]   = column;
          J_vals_buffer_[nnz_++] = value;
        }
      return entries == 0 ? 0 : this->constructCoo();
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

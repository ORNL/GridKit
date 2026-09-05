#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/EMT/PhasorInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase impedance load
     *
     * System sizes:
     * - Number of equations = 3
     * - Number of independent variables = 3
     */
    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ()
      : LoadZ(ModelDataT{})
    {
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      equation_size_ = 3;
      supplied_fit_  = data.Z.has_value();
      if (supplied_fit_)
      {
        if (data.Z->rows != 3 || data.Z->cols != 3)
        {
          throw std::invalid_argument("LoadZ: Z must be a three-phase operator");
        }
        z_.emplace(*data.Z, ONE<RealT>);
        derivative_columns_independent_ = data.Z->E.hasShape(3, 3)
                                          && independentDerivativeColumns<RealT>(data.Z->E);
      }
      else
      {
        VectorFitData<RealT, IdxT> impedance;
        impedance.D = R_;
        impedance.E = L_;
        z_.emplace(impedance, ONE<RealT>);
        derivative_columns_independent_ = independentDerivativeColumns<RealT>(L_);
      }
      this->addOperator(&*z_);
      size_ = equation_size_ + z_->size();
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::~LoadZ()
    {
    }

    /**
     * @brief Read model parameters from the data object
     */
    template <typename scalar_type, typename index_type>
    void LoadZ<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::N))
      {
        n_phases_ = std::get<IdxT>(data.parameters.at(Parameter::N));
      }

      if (data.parameters.contains(Parameter::R))
      {
        R_ = std::get<ABCMatrix<RealT>>(data.parameters.at(Parameter::R));
      }

      if (data.parameters.contains(Parameter::L))
      {
        L_ = std::get<ABCMatrix<RealT>>(data.parameters.at(Parameter::L));
      }
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local storage and registers coupling signals.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Bind the current port and wire the rational impedance before its
      // allocation, so index assignment can route into it
      this->bindPort(i_port_, 0);
      z_->attachInput(&i_port_);
      z_->attachOutput(&i_port_);
      const int status = this->allocateOperators();
      if (status != 0)
      {
        return status;
      }

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize coupling data
      this->allocateExternalVectors(static_cast<IdxT>(LoadZExternalVariables::MAXIMUM), 3);
      signals_.registerExternalVariableSignals(*this);
      this->setExternalResidualSignal(0, signals_.template getAttachedSignal<LoadZExternalVariables::VA>());
      this->setExternalResidualSignal(1, signals_.template getAttachedSignal<LoadZExternalVariables::VB>());
      this->setExternalResidualSignal(2, signals_.template getAttachedSignal<LoadZExternalVariables::VC>());

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     *
     * @return Number of model configuration errors found
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::verify() const
    {
      int error_count = z_->verify();
      if (!derivative_columns_independent_)
      {
        Log::error() << "LoadZ: nonzero impedance derivative columns must be independent\n";
        ++error_count;
      }
      if (n_phases_ != IdxT{3})
      {
        Log::error() << "LoadZ: only three phases are supported\n";
        ++error_count;
      }

      if (!signals_.template isAttached<LoadZExternalVariables::VA>()
          || !signals_.template isAttached<LoadZExternalVariables::VB>()
          || !signals_.template isAttached<LoadZExternalVariables::VC>())
      {
        Log::error() << "LoadZ: the bus voltage port is not attached\n";
        ++error_count;
      }

      if (supplied_fit_ && (hasResistance() || hasInductance()))
      {
        Log::error() << "LoadZ: the rational impedance excludes the "
                        "resistance and inductance matrices\n";
        ++error_count;
      }

      return error_count;
    }

    /**
     * Initialization of the load model
     *
     * A purely resistive load solves its algebraic current from the bus
     * voltage. An inductive load starts de-energized and the integrator
     * initialization establishes consistent conditions.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::initialize()
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = 0.0;
      y[1]  = 0.0;
      y[2]  = 0.0;
      yp[0] = 0.0;
      yp[1] = 0.0;
      yp[2] = 0.0;

      if (!hasInductance() && !supplied_fit_)
      {
        this->gatherExternalVariables();
        ABCVector<RealT> voltage{};
        ABCVector<RealT> current{};
        ABCVector<RealT> derivative{};
        for (size_t n = 0; n < 3; ++n)
        {
          voltage[n] = -static_cast<RealT>(y_ext_[n]);
        }
        if (solvePhasorSystem(R_, ABCMatrix<RealT>{}, ONE<RealT>, voltage, ABCVector<RealT>{}, current, derivative) == 0)
        {
          for (size_t n = 0; n < 3; ++n)
          {
            y[n] = current[n];
          }
        }
        // Singular resistance leaves the zero-current guess for the
        // assembled circuit initialization, including an ideal short.
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return this->initializeOperators();
    }

    /**
     * @brief Solve Z(j omega) I = -V and initialize the impedance memory.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::initializeSteadyState(RealT omega)
    {
      if (!std::isfinite(omega) || omega <= ZERO<RealT>)
      {
        return 1;
      }
      this->gatherExternalVariables();
      ABCVector<RealT> v{};
      ABCVector<RealT> v_dot{};
      ABCVector<RealT> current{};
      ABCVector<RealT> current_dot{};
      for (size_t n = 0; n < 3; ++n)
      {
        v[n]     = -static_cast<RealT>(y_ext_[n]);
        v_dot[n] = -static_cast<RealT>(yp_ext_[n]);
      }
      ABCMatrix<RealT> Z_re{};
      ABCMatrix<RealT> Z_im{};
      try
      {
        z_->transfer(omega, Z_re, Z_im);
      }
      catch (const std::domain_error&)
      {
        return 1;
      }
      if (solvePhasorSystem(Z_re, Z_im, omega, v, v_dot, current, current_dot) != 0)
      {
        return 1;
      }
      auto* y  = y_.getData();
      auto* yp = yp_.getData();
      for (size_t n = 0; n < 3; ++n)
      {
        y[n]  = current[n];
        yp[n] = current_dot[n];
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return z_->initializeSteadyState(omega, current, current_dot);
    }

    /**
     * @brief Whether any resistance entry is nonzero.
     */
    template <typename scalar_type, typename index_type>
    bool LoadZ<scalar_type, index_type>::hasResistance() const
    {
      bool nonzero = false;
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          if (R_[n][k] != 0.0)
          {
            nonzero = true;
          }
        }
      }
      return nonzero;
    }

    /**
     * @brief Whether any inductance entry is nonzero.
     *
     * Used to retain the de-energized startup of legacy inductive loads.
     */
    template <typename scalar_type, typename index_type>
    bool LoadZ<scalar_type, index_type>::hasInductance() const
    {
      bool nonzero = false;
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          if (L_[n][k] != 0.0)
          {
            nonzero = true;
          }
        }
      }
      return nonzero;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::tagDifferentiable()
    {
      for (size_t n = 0; n < 3; ++n)
      {
        tag_[n] = z_->hasInputDerivative(static_cast<IdxT>(n));
      }
      return this->tagDifferentiableOperators();
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
    int LoadZ<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return this->setAbsoluteToleranceOperators(rel_tol);
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZ<scalar_type, index_type>::evaluateInternalResidual(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      // The impedance output accumulates into these branch rows.
      f[0] = y_ext[0];
      f[1] = y_ext[1];
      f[2] = y_ext[2];

      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZ<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT ia = y[0];
      const ScalarT ib = y[1];
      const ScalarT ic = y[2];

      f_ext[0] = ia;
      f_ext[1] = ib;
      f_ext[2] = ic;

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();

      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      evaluateInternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f);
      this->evaluateOperatorInternalResiduals();
      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief External residual contributions to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateExternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f_ext_.data());
      this->scatterExternalResidual();
      this->evaluateOperatorExternalResiduals();

      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
    }

    /**
     * @brief Exact local Jacobian with signal chain rules and operator contributions.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateJacobian()
    {
      this->gatherExternalVariables();
      const auto&                              external = this->externalVariableSignals();
      std::vector<typename SignalT::GradientT> gradients(external.size());
      size_t                                   capacity = 3;
      for (size_t k = 0; k < external.size(); ++k)
      {
        external[k]->appendGradient(gradients[k]);
        capacity += gradients[k].size();
      }
      capacity += static_cast<size_t>(z_->jacobianCapacity());
      if (this->hasComputedInputs() || capacity > jacobian_capacity_)
      {
        this->resetJacobianStructure();
      }
      if (capacity > jacobian_capacity_)
      {
        delete[] J_rows_buffer_;
        delete[] J_cols_buffer_;
        delete[] J_vals_buffer_;
        J_rows_buffer_     = new IdxT[capacity];
        J_cols_buffer_     = new IdxT[capacity];
        J_vals_buffer_     = new RealT[capacity];
        jacobian_capacity_ = capacity;
      }
      nnz_        = 0;
      auto append = [this](IdxT row, IdxT column, RealT value)
      {
        J_rows_buffer_[nnz_]   = row;
        J_cols_buffer_[nnz_]   = column;
        J_vals_buffer_[nnz_++] = value;
      };
      for (size_t n = 0; n < 3; ++n)
      {
        append(residual_indices_ext_[n], variable_indices_[n], ONE<RealT>);
        for (const auto& [column, value] : gradients[n])
        {
          append(residual_indices_[n], column, value);
        }
      }
      const int status = this->evaluateOperatorJacobians();
      if (status != 0)
        return status;
      this->appendOperatorJacobians();
      return this->constructCoo();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* LoadZ<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void LoadZ<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::ia, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::ib, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::ic, [this]
                    { return y_.getData()[2]; });
    }

  } // namespace EMT
} // namespace GridKit

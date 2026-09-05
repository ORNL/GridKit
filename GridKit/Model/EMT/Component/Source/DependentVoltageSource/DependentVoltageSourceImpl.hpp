#pragma once

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSource.hpp>
#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSourceData.hpp>
#include <GridKit/Model/EMT/PhasorInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase dependent voltage source
     *
     * System sizes:
     * - Number of equations = 3
     * - Number of independent variables = 3
     */
    template <typename scalar_type, typename index_type>
    DependentVoltageSource<scalar_type, index_type>::DependentVoltageSource()
      : DependentVoltageSource(ModelDataT{})
    {
    }

    template <typename scalar_type, typename index_type>
    DependentVoltageSource<scalar_type, index_type>::DependentVoltageSource(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      equation_size_ = 3;
      if (data.Y.has_value())
      {
        if (data.Y->rows != 3 || data.Y->cols != 3)
        {
          throw std::invalid_argument("DependentVoltageSource: Y must be a three-phase operator");
        }
        yfit_.emplace(*data.Y, ONE<RealT>);
        this->addOperator(&*yfit_);
        rl_on_          = ZERO<RealT>;
        fit_on_         = ONE<RealT>;
        fit_ey_nonzero_ = yfit_->hasFeedthroughDerivative();
      }
      else
      {
        VectorFitData<RealT, IdxT> impedance;
        impedance.D = Rs_;
        impedance.E = Ls_;
        zfit_.emplace(impedance, ONE<RealT>);
        this->addOperator(&*zfit_);
      }
      size_ = equation_size_ + (yfit_.has_value() ? yfit_->size() : zfit_->size());
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    DependentVoltageSource<scalar_type, index_type>::~DependentVoltageSource()
    {
    }

    /**
     * @brief Read model parameters from the data object
     */
    template <typename scalar_type, typename index_type>
    void DependentVoltageSource<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::N))
      {
        n_phases_ = std::get<IdxT>(data.parameters.at(Parameter::N));
      }

      if (data.parameters.contains(Parameter::Rs))
      {
        Rs_ = std::get<ABCMatrix<RealT>>(data.parameters.at(Parameter::Rs));
      }

      if (data.parameters.contains(Parameter::Ls))
      {
        Ls_ = std::get<ABCMatrix<RealT>>(data.parameters.at(Parameter::Ls));
      }
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local storage and registers coupling signals.
     */
    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Bind the branch voltage port and wire the rational admittance before
      // its allocation, so index assignment can route into it
      this->bindPort(u_port_, 0);
      if (yfit_.has_value())
      {
        yfit_->attachInput(&u_port_);
        yfit_->attachOutput(signals_.template getAttachedSignal<DependentVoltageSourceExternalVariables::VA>(),
                            signals_.template getAttachedSignal<DependentVoltageSourceExternalVariables::VB>(),
                            signals_.template getAttachedSignal<DependentVoltageSourceExternalVariables::VC>());
      }
      if (zfit_.has_value())
      {
        zfit_->attachInput(&u_port_);
        zfit_->attachOutput(&u_port_);
      }
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
      this->allocateExternalVectors(static_cast<IdxT>(DependentVoltageSourceExternalVariables::MAXIMUM), 3);
      signals_.registerExternalVariableSignals(*this);
      this->setExternalResidualSignal(0, signals_.template getAttachedSignal<DependentVoltageSourceExternalVariables::VA>());
      this->setExternalResidualSignal(1, signals_.template getAttachedSignal<DependentVoltageSourceExternalVariables::VB>());
      this->setExternalResidualSignal(2, signals_.template getAttachedSignal<DependentVoltageSourceExternalVariables::VC>());

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     *
     * @return Number of model configuration errors found
     */
    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::verify() const
    {
      int error_count = 0;

      if (zfit_.has_value())
      {
        error_count += zfit_->verify();
        if (!independentDerivativeColumns<RealT>(Ls_))
        {
          Log::error() << "DependentVoltageSource: nonzero series-inductance columns must be independent\n";
          ++error_count;
        }
      }

      if (n_phases_ != IdxT{3})
      {
        Log::error() << "DependentVoltageSource: only three phases are supported\n";
        ++error_count;
      }

      if (!signals_.template isAttached<DependentVoltageSourceExternalVariables::VA>()
          || !signals_.template isAttached<DependentVoltageSourceExternalVariables::VB>()
          || !signals_.template isAttached<DependentVoltageSourceExternalVariables::VC>())
      {
        Log::error() << "DependentVoltageSource: the bus voltage port is not attached\n";
        ++error_count;
      }

      if (!signals_.template isAttached<DependentVoltageSourceExternalVariables::EA>()
          || !signals_.template isAttached<DependentVoltageSourceExternalVariables::EB>()
          || !signals_.template isAttached<DependentVoltageSourceExternalVariables::EC>())
      {
        Log::error() << "DependentVoltageSource: the source voltage input is not attached\n";
        ++error_count;
      }
      else if (!signals_.template isLinked<DependentVoltageSourceExternalVariables::EA>()
               || !signals_.template isLinked<DependentVoltageSourceExternalVariables::EB>()
               || !signals_.template isLinked<DependentVoltageSourceExternalVariables::EC>())
      {
        Log::error() << "DependentVoltageSource: the source voltage input has no linked source\n";
        ++error_count;
      }

      if (yfit_.has_value())
      {
        error_count += yfit_->verify();

        bool matrices_nonzero = false;
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            if (Rs_[n][k] != 0.0 || Ls_[n][k] != 0.0)
            {
              matrices_nonzero = true;
            }
          }
        }
        if (matrices_nonzero)
        {
          Log::error() << "DependentVoltageSource: the rational admittance excludes "
                          "the series matrices\n";
          ++error_count;
        }

        if (fit_ey_nonzero_)
        {
          Log::error() << "DependentVoltageSource: the rational admittance linear "
                          "coefficient must be zero, because the branch "
                          "voltage is algebraic\n";
          ++error_count;
        }
      }

      return error_count;
    }

    /**
     * Initialization of the dependent voltage source model
     *
     */
    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::initialize()
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0] = 0.0;
      y[1] = 0.0;
      y[2] = 0.0;

      yp[0] = 0.0;
      yp[1] = 0.0;
      yp[2] = 0.0;

      if (yfit_.has_value())
      {
        this->gatherExternalVariables();
        for (size_t n = 0; n < 3; ++n)
        {
          y[n]  = y_ext_[3 + n] - y_ext_[n];
          yp[n] = yp_ext_[3 + n] - yp_ext_[n];
        }
      }
      const int status = this->initializeOperators();

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return status;
    }

    /**
     * @brief Initialize branch samples and rational memory in sinusoidal steady state.
     */
    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::initializeSteadyState(RealT omega)
    {
      if (!std::isfinite(omega) || omega <= ZERO<RealT>)
      {
        return 1;
      }
      this->gatherExternalVariables();
      auto*            y  = y_.getData();
      auto*            yp = yp_.getData();
      ABCVector<RealT> u{};
      ABCVector<RealT> u_dot{};
      for (size_t n = 0; n < 3; ++n)
      {
        u[n]     = static_cast<RealT>(y_ext_[3 + n] - y_ext_[n]);
        u_dot[n] = static_cast<RealT>(yp_ext_[3 + n] - yp_ext_[n]);
      }
      ABCVector<RealT> branch     = u;
      ABCVector<RealT> branch_dot = u_dot;
      if (zfit_.has_value())
      {
        ABCMatrix<RealT> Z_re{};
        ABCMatrix<RealT> Z_im{};
        zfit_->transfer(omega, Z_re, Z_im);
        if (solvePhasorSystem(Z_re, Z_im, omega, u, u_dot, branch, branch_dot) != 0)
        {
          return 1;
        }
      }
      for (size_t n = 0; n < 3; ++n)
      {
        y[n]  = branch[n];
        yp[n] = branch_dot[n];
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return yfit_.has_value() ? yfit_->initializeSteadyState(omega, branch, branch_dot)
                               : zfit_->initializeSteadyState(omega, branch, branch_dot);
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::tagDifferentiable()
    {
      for (size_t n = 0; n < 3; ++n)
      {
        tag_[n] = zfit_.has_value() && zfit_->hasInputDerivative(static_cast<IdxT>(n));
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
    int DependentVoltageSource<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return this->setAbsoluteToleranceOperators(rel_tol);
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int DependentVoltageSource<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      /* Read variables */
      const ScalarT ia = y[0];
      const ScalarT ib = y[1];
      const ScalarT ic = y[2];

      // Set coupling variable aliases
      const ScalarT va = y_ext[0];
      const ScalarT vb = y_ext[1];
      const ScalarT vc = y_ext[2];
      const ScalarT ea = y_ext[3];
      const ScalarT eb = y_ext[4];
      const ScalarT ec = y_ext[5];

      // The admittance input is the branch voltage. The legacy impedance
      // contributes its voltage drop through its VectorFit output.
      f[0] = fit_on_ * ia + va - ea;
      f[1] = fit_on_ * ib + vb - eb;
      f[2] = fit_on_ * ic + vc - ec;

      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int DependentVoltageSource<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT ia = y[0];
      const ScalarT ib = y[1];
      const ScalarT ic = y[2];

      f_ext[0] = rl_on_ * ia;
      f_ext[1] = rl_on_ * ib;
      f_ext[2] = rl_on_ * ic;

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::evaluateInternalResidual()
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
    int DependentVoltageSource<scalar_type, index_type>::evaluateExternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f_ext_.data());
      this->scatterExternalResidual();
      this->evaluateOperatorExternalResiduals();

      return 0;
    }

    /**
     * @brief Residual contribution of the source is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
    }

    /**
     * @brief Exact local Jacobian with signal chain rules and operator contributions.
     */
    template <typename scalar_type, typename index_type>
    int DependentVoltageSource<scalar_type, index_type>::evaluateJacobian()
    {
      this->gatherExternalVariables();
      const auto&                              external = this->externalVariableSignals();
      std::vector<typename SignalT::GradientT> gradients(external.size());
      size_t                                   capacity = 6;
      for (size_t k = 0; k < external.size(); ++k)
      {
        external[k]->appendGradient(gradients[k]);
        capacity += gradients[k].size();
      }
      if (yfit_.has_value())
        capacity += static_cast<size_t>(yfit_->jacobianCapacity());
      if (zfit_.has_value())
        capacity += static_cast<size_t>(zfit_->jacobianCapacity());
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
        append(residual_indices_[n], variable_indices_[n], fit_on_);
        append(residual_indices_ext_[n], variable_indices_[n], rl_on_);
        for (const auto& [column, value] : gradients[n])
        {
          append(residual_indices_[n], column, value);
        }
        for (const auto& [column, value] : gradients[3 + n])
        {
          append(residual_indices_[n], column, -value);
        }
      }
      const int status = this->evaluateOperatorJacobians();
      if (status != 0)
        return status;
      this->appendOperatorJacobians();
      return this->constructCoo();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* DependentVoltageSource<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void DependentVoltageSource<scalar_type, index_type>::initializeMonitor()
    {
      using Variable  = typename ModelDataT::MonitorableVariables;
      using Variables = DependentVoltageSourceExternalVariables;

      monitor_->set(Variable::ea, [this]
                    { return signals_.template readExternalVariable<Variables::EA>(); });
      monitor_->set(Variable::eb, [this]
                    { return signals_.template readExternalVariable<Variables::EB>(); });
      monitor_->set(Variable::ec, [this]
                    { return signals_.template readExternalVariable<Variables::EC>(); });
      monitor_->set(Variable::ia, [this]
                    { return yfit_.has_value() ? yfit_->output(0) : y_.getData()[0]; });
      monitor_->set(Variable::ib, [this]
                    { return yfit_.has_value() ? yfit_->output(1) : y_.getData()[1]; });
      monitor_->set(Variable::ic, [this]
                    { return yfit_.has_value() ? yfit_->output(2) : y_.getData()[2]; });
    }

  } // namespace EMT
} // namespace GridKit

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>

#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumped.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedData.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase lumped line
     *
     * System sizes:
     * - Number of equations = 3
     * - Number of independent variables = 3
     */
    template <typename scalar_type, typename index_type>
    LineLumped<scalar_type, index_type>::LineLumped()
    {
      equation_size_ = size_ = 3;
      setDerivedParams();
      initializePorts();
    }

    template <typename scalar_type, typename index_type>
    LineLumped<scalar_type, index_type>::LineLumped(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      equation_size_ = size_ = 3;
      setDerivedParams();
      if (data.Zp.has_value() && data.Yp.has_value())
      {
        z_.emplace(*data.Zp, dx_);
        this->addOperator(&*z_);
        size_  = equation_size_ + z_->size();
        rl_on_ = ZERO<RealT>;

        // The series linear coefficient must be nonsingular so the series
        // current stays differential
        const auto& E   = data.Zp->E;
        const RealT det = E[0][0] * (E[1][1] * E[2][2] - E[1][2] * E[2][1])
                          - E[0][1] * (E[1][0] * E[2][2] - E[1][2] * E[2][0])
                          + E[0][2] * (E[1][0] * E[2][1] - E[1][1] * E[2][0]);
        fit_ez_singular_ = det == 0.0;
      }
      initializePorts();
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    LineLumped<scalar_type, index_type>::~LineLumped()
    {
    }

    /**
     * @brief Read model parameters from the data object
     */
    template <typename scalar_type, typename index_type>
    void LineLumped<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::conductors))
      {
        conductors_ = std::get<ABCVector<IdxT>>(data.parameters.at(Parameter::conductors));
      }

      if (data.parameters.contains(Parameter::dx))
      {
        dx_ = std::get<RealT>(data.parameters.at(Parameter::dx));
      }

      if (data.parameters.contains(Parameter::Rp))
      {
        Rp_ = std::get<ABCMatrix<RealT>>(data.parameters.at(Parameter::Rp));
      }

      if (data.parameters.contains(Parameter::Lp))
      {
        Lp_ = std::get<ABCMatrix<RealT>>(data.parameters.at(Parameter::Lp));
      }

      if (data.parameters.contains(Parameter::Gp))
      {
        Gp_ = std::get<ABCMatrix<RealT>>(data.parameters.at(Parameter::Gp));
      }

      if (data.parameters.contains(Parameter::Cp))
      {
        Cp_ = std::get<ABCMatrix<RealT>>(data.parameters.at(Parameter::Cp));
      }
    }

    /**
     * @brief Derived parameters
     *
     * The per-unit-length matrices are scaled by the segment length.
     */
    template <typename scalar_type, typename index_type>
    void LineLumped<scalar_type, index_type>::setDerivedParams()
    {
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          R_[n][k] = dx_ * Rp_[n][k];
          L_[n][k] = dx_ * Lp_[n][k];
        }
      }
    }

    template <typename scalar_type, typename index_type>
    void LineLumped<scalar_type, index_type>::initializePorts()
    {
      for (size_t p = 0; p < 3; ++p)
        i21_port_[p].setComputed(
            [this, p]
            { return -i12_port_[p].read(); },
            [this, p](typename SignalT::GradientT& gradient, RealT scale)
            { i12_port_[p].appendGradient(gradient, -scale); });
    }

    template <typename scalar_type, typename index_type>
    void LineLumped<scalar_type, index_type>::attachTerminal(size_t end, PhaseSignals voltage)
    {
      if (allocated_)
        throw std::logic_error("LineLumped terminals cannot change after allocation");
      if (end > 1)
        throw std::out_of_range("Invalid LineLumped terminal");
      for (size_t p = 0; p < 3; ++p)
        if (!voltage[p])
          throw std::invalid_argument("LineLumped requires terminal voltage inputs");
      for (size_t p = 0; p < 3; ++p)
        signals_.attachSignal(static_cast<LineLumpedExternalVariables>(3 * end + p), voltage[p]);
    }

    template <typename scalar_type, typename index_type>
    typename LineLumped<scalar_type, index_type>::SignalT& LineLumped<scalar_type, index_type>::outputSignal(Outputs output)
    {
      const size_t index = static_cast<size_t>(output);
      if (index < 3)
        return i12_port_[index];
      return i21_port_.at(index - 3);
    }

    template <typename scalar_type, typename index_type>
    void LineLumped<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      if (static_cast<size_t>(output) < 3)
        signals_.assignSignal(static_cast<LineLumpedInternalVariables>(output), signal);
      else
      {
        auto* value = &outputSignal(output);
        signal->claimProducer();
        signal->setComputed([value]
                            { return value->read(); },
                            [value](typename SignalT::GradientT& gradient, RealT scale)
                            { value->appendGradient(gradient, scale); });
      }
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local storage and registers coupling signals.
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      for (IdxT phase = 0; phase < 3; ++phase)
        this->bindSignal(i12_port_[static_cast<size_t>(phase)], phase);
      if (z_.has_value())
      {
        z_->attachInput(&i12_port_[0], &i12_port_[1], &i12_port_[2]);
        z_->attachOutput(&i12_port_[0], &i12_port_[1], &i12_port_[2]);
        const int status = this->allocateOperators();
        if (status != 0)
          return status;
      }

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      this->allocateExternalVectors(static_cast<IdxT>(LineLumpedExternalVariables::MAXIMUM), 0);
      signals_.registerExternalVariableSignals(*this);
      signals_.bindInternalVariableSignals(*this);
      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     *
     * @return Number of model configuration errors found
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::verify() const
    {
      int error_count = 0;

      if (!signals_.template isAttached<LineLumpedExternalVariables::V1A>()
          || !signals_.template isAttached<LineLumpedExternalVariables::V1B>()
          || !signals_.template isAttached<LineLumpedExternalVariables::V1C>()
          || !signals_.template isAttached<LineLumpedExternalVariables::V2A>()
          || !signals_.template isAttached<LineLumpedExternalVariables::V2B>()
          || !signals_.template isAttached<LineLumpedExternalVariables::V2C>())
      {
        Log::error() << "LineLumped: a terminal voltage port is not attached\n";
        ++error_count;
      }

      if (dx_ <= 0.0)
      {
        Log::error() << "LineLumped: the segment length must be positive\n";
        ++error_count;
      }

      if (conductors_[0] != 1 || conductors_[1] != 2 || conductors_[2] != 3)
      {
        Log::error() << "LineLumped: the conductor phase-index list must be [1, 2, 3]\n";
        ++error_count;
      }

      if (z_.has_value())
      {
        error_count           += z_->verify();
        bool matrices_nonzero  = false;
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            if (Rp_[n][k] != 0.0 || Lp_[n][k] != 0.0 || Gp_[n][k] != 0.0 || Cp_[n][k] != 0.0)
            {
              matrices_nonzero = true;
            }
          }
        }
        if (matrices_nonzero)
        {
          Log::error() << "LineLumped: the rational operators exclude the "
                          "per-unit-length matrices\n";
          ++error_count;
        }

        if (fit_ez_singular_)
        {
          Log::error() << "LineLumped: the series linear coefficient must be "
                          "nonsingular\n";
          ++error_count;
        }
      }

      return error_count;
    }

    /**
     * Initialization of the lumped line model
     *
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      for (const auto& [output, value] : outputs)
        if (static_cast<size_t>(output) >= 3)
          throw std::invalid_argument("LineLumped initial outputs must be series currents");
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      for (IdxT j = 0; j < this->equationSize(); ++j)
      {
        y[j]  = 0.0;
        yp[j] = 0.0;
      }

      for (const auto& [output, value] : outputs)
      {
        y[static_cast<size_t>(output)] = static_cast<ScalarT>(value);
      }

      if (z_.has_value())
      {
        const int status = z_->initialize();
        if (status != 0)
          return status;
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

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
    int LineLumped<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return this->setAbsoluteToleranceOperators(rel_tol);
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LineLumped<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      /* Read variables */
      const ScalarT i12a = y[0];
      const ScalarT i12b = y[1];
      const ScalarT i12c = y[2];

      /* Read derivatives */
      const ScalarT i12a_dot = yp[0];
      const ScalarT i12b_dot = yp[1];
      const ScalarT i12c_dot = yp[2];

      // Set coupling variable aliases
      const ScalarT v1a = y_ext[0];
      const ScalarT v1b = y_ext[1];
      const ScalarT v1c = y_ext[2];
      const ScalarT v2a = y_ext[3];
      const ScalarT v2b = y_ext[4];
      const ScalarT v2c = y_ext[5];

      /* 3 series branch equations; the rational series terms accumulate
         through the operator when the matrix mask is off */
      f[0] = rl_on_
                 * (R_[0][0] * i12a + R_[0][1] * i12b + R_[0][2] * i12c
                    + L_[0][0] * i12a_dot + L_[0][1] * i12b_dot + L_[0][2] * i12c_dot)
             + v2a - v1a;
      f[1] = rl_on_
                 * (R_[1][0] * i12a + R_[1][1] * i12b + R_[1][2] * i12c
                    + L_[1][0] * i12a_dot + L_[1][1] * i12b_dot + L_[1][2] * i12c_dot)
             + v2b - v1b;
      f[2] = rl_on_
                 * (R_[2][0] * i12a + R_[2][1] * i12b + R_[2][2] * i12c
                    + L_[2][0] * i12a_dot + L_[2][1] * i12b_dot + L_[2][2] * i12c_dot)
             + v2c - v1c;

      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LineLumped<scalar_type, index_type>::evaluateExternalResidual(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        [[maybe_unused]] ScalarT*       f_ext)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::evaluateInternalResidual()
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
     * @brief Add the embedded series-impedance contribution.
     *
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::evaluateExternalResidual()
    {
      return this->evaluateOperatorExternalResiduals();
    }

    /**
     * @brief Evaluate the series-current equation and its embedded impedance.
     *
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* LineLumped<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void LineLumped<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::i12a, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::i12b, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::i12c, [this]
                    { return y_.getData()[2]; });
    }

  } // namespace EMT
} // namespace GridKit

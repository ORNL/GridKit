#pragma once

#include <algorithm>
#include <iostream>

#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZData.hpp>
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
    {
      size_ = 3;
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      size_ = 3;
      if (data.Z.has_value())
      {
        z_.emplace(*data.Z, ONE<RealT>);
        equation_size_ = 3;
        this->addOperator(&*z_);
        size_  = equation_size_ + z_->size();
        rl_on_ = ZERO<RealT>;

        // The specification admits a zero or nonsingular linear coefficient
        const auto& E   = data.Z->E;
        const RealT det = E[0][0] * (E[1][1] * E[2][2] - E[1][2] * E[2][1])
                          - E[0][1] * (E[1][0] * E[2][2] - E[1][2] * E[2][0])
                          + E[0][2] * (E[1][0] * E[2][1] - E[1][1] * E[2][0]);
        fit_e_singular_ = z_->hasFeedthroughDerivative() && det == 0.0;
      }
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
      if (z_.has_value())
      {
        z_->attachInput(&i_port_);
        z_->attachOutput(&i_port_);
        this->allocateOperators();
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
      int error_count = 0;

      if (!signals_.template isAttached<LoadZExternalVariables::VA>()
          || !signals_.template isAttached<LoadZExternalVariables::VB>()
          || !signals_.template isAttached<LoadZExternalVariables::VC>())
      {
        Log::error() << "LoadZ: the bus voltage port is not attached\n";
        ++error_count;
      }

      if (z_.has_value())
      {
        error_count += z_->verify();

        if (hasResistance() || hasInductance())
        {
          Log::error() << "LoadZ: the rational impedance excludes the "
                          "resistance and inductance matrices\n";
          ++error_count;
        }

        if (fit_e_singular_)
        {
          Log::error() << "LoadZ: the rational impedance linear coefficient "
                          "must be zero or nonsingular\n";
          ++error_count;
        }
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

      if (z_.has_value())
      {
        this->initializeOperators();
      }

      if (!hasInductance() && !z_.has_value())
      {
        const ScalarT va = signals_.template readExternalVariable<LoadZExternalVariables::VA>();
        const ScalarT vb = signals_.template readExternalVariable<LoadZExternalVariables::VB>();
        const ScalarT vc = signals_.template readExternalVariable<LoadZExternalVariables::VC>();

        // Solve R i = -v by the cofactor inverse
        const RealT det = R_[0][0] * (R_[1][1] * R_[2][2] - R_[1][2] * R_[2][1])
                          - R_[0][1] * (R_[1][0] * R_[2][2] - R_[1][2] * R_[2][0])
                          + R_[0][2] * (R_[1][0] * R_[2][1] - R_[1][1] * R_[2][0]);
        if (det != 0.0)
        {
          y[0] = -((R_[1][1] * R_[2][2] - R_[1][2] * R_[2][1]) * va
                   + (R_[0][2] * R_[2][1] - R_[0][1] * R_[2][2]) * vb
                   + (R_[0][1] * R_[1][2] - R_[0][2] * R_[1][1]) * vc)
                 / det;
          y[1] = -((R_[1][2] * R_[2][0] - R_[1][0] * R_[2][2]) * va
                   + (R_[0][0] * R_[2][2] - R_[0][2] * R_[2][0]) * vb
                   + (R_[0][2] * R_[1][0] - R_[0][0] * R_[1][2]) * vc)
                 / det;
          y[2] = -((R_[1][0] * R_[2][1] - R_[1][1] * R_[2][0]) * va
                   + (R_[0][1] * R_[2][0] - R_[0][0] * R_[2][1]) * vb
                   + (R_[0][0] * R_[1][1] - R_[0][1] * R_[1][0]) * vc)
                 / det;
        }
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
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
     * The injected current is differential when the inductance matrix is
     * nonsingular and algebraic when it is zero.
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
      bool differential = hasInductance();
      if (z_.has_value())
      {
        differential = z_->hasFeedthroughDerivative();
      }

      tag_[0] = differential;
      tag_[1] = differential;
      tag_[2] = differential;

      if (z_.has_value())
      {
        this->tagDifferentiableOperators();
      }

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
    int LoadZ<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZ<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      /* Read variables */
      const ScalarT ia = y[0];
      const ScalarT ib = y[1];
      const ScalarT ic = y[2];

      /* Read derivatives */
      const ScalarT ia_dot = yp[0];
      const ScalarT ib_dot = yp[1];
      const ScalarT ic_dot = yp[2];

      // Set coupling variable aliases
      const ScalarT va = y_ext[0];
      const ScalarT vb = y_ext[1];
      const ScalarT vc = y_ext[2];

      /* 3 load branch equations; the rational impedance terms accumulate
         through the operator when the resistance and inductance mask is off */
      f[0] = rl_on_
                 * (R_[0][0] * ia + R_[0][1] * ib + R_[0][2] * ic
                    + L_[0][0] * ia_dot + L_[0][1] * ib_dot + L_[0][2] * ic_dot)
             + va;
      f[1] = rl_on_
                 * (R_[1][0] * ia + R_[1][1] * ib + R_[1][2] * ic
                    + L_[1][0] * ia_dot + L_[1][1] * ib_dot + L_[1][2] * ic_dot)
             + vb;
      f[2] = rl_on_
                 * (R_[2][0] * ia + R_[2][1] * ib + R_[2][2] * ic
                    + L_[2][0] * ia_dot + L_[2][1] * ib_dot + L_[2][2] * ic_dot)
             + vc;

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

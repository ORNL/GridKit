#pragma once

#include <algorithm>
#include <iostream>

#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumped.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase lumped line
     *
     * System sizes:
     * - Number of equations = 9
     * - Number of independent variables = 9
     */
    template <typename scalar_type, typename index_type>
    LineLumped<scalar_type, index_type>::LineLumped()
    {
      size_ = 9;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LineLumped<scalar_type, index_type>::LineLumped(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      size_ = 9;
      setDerivedParams();
      if (data.Zp.has_value() && data.Yp.has_value())
      {
        own_size_ = 9;
        z_.emplace(*data.Zp, dx_);
        this->registerSubmodel(&*z_);
        y1_.emplace(*data.Yp, dx_);
        this->registerSubmodel(&*y1_);
        y2_.emplace(*data.Yp, dx_);
        this->registerSubmodel(&*y2_);
        size_  = own_size_ + z_->size() + y1_->size() + y2_->size();
        rl_on_ = ZERO<RealT>;

        // The series linear coefficient must be nonsingular so the series
        // current stays differential
        const auto& E   = data.Zp->E;
        const RealT det = E[0][0] * (E[1][1] * E[2][2] - E[1][2] * E[2][1])
                          - E[0][1] * (E[1][0] * E[2][2] - E[1][2] * E[2][0])
                          + E[0][2] * (E[1][0] * E[2][1] - E[1][1] * E[2][0]);
        fit_ez_singular_ = det == 0.0;
      }
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
          G_[n][k] = dx_ * Gp_[n][k];
          C_[n][k] = dx_ * Cp_[n][k];
        }
      }
    }

    /**
     * @brief Whether any shunt capacitance entry is nonzero.
     *
     * A nonzero shunt capacitance makes the shunt rows read the terminal
     * voltage derivatives.
     */
    template <typename scalar_type, typename index_type>
    bool LineLumped<scalar_type, index_type>::hasShuntCapacitance() const
    {
      bool nonzero = false;
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          if (C_[n][k] != 0.0)
          {
            nonzero = true;
          }
        }
      }
      return nonzero;
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
     * @brief allocate method resizes local storage and registers coupling nodes.
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

      // Bind the series current and shunt-row ports and wire the rational
      // submodels before their allocation, so index assignment can route
      // into them
      this->bindPort(i12_port_, 0);
      this->bindPort(sh1_rows_port_, 3);
      this->bindPort(sh2_rows_port_, 6);
      if (z_.has_value())
      {
        z_->attachInput(&i12_port_);
        z_->attachOutput(&i12_port_);
        y1_->attachInput(signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V1A>(),
                         signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V1B>(),
                         signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V1C>());
        y1_->attachOutput(&sh1_rows_port_);
        y2_->attachInput(signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V2A>(),
                         signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V2B>(),
                         signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V2C>());
        y2_->attachOutput(&sh2_rows_port_);
        this->allocateSubmodels();
      }

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize coupling data
      this->allocateExternalVectors(static_cast<IdxT>(LineLumpedExternalVariables::MAXIMUM), 6);
      signals_.registerExternalVariableNodes(*this);
      this->setExternalResidualNode(0, signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V1A>());
      this->setExternalResidualNode(1, signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V1B>());
      this->setExternalResidualNode(2, signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V1C>());
      this->setExternalResidualNode(3, signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V2A>());
      this->setExternalResidualNode(4, signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V2B>());
      this->setExternalResidualNode(5, signals_.template getAttachedSignalNode<LineLumpedExternalVariables::V2C>());

      // The shunt rows read the terminal voltage derivatives, so the
      // connected bus voltages become differential.
      if (hasShuntCapacitance())
      {
        signals_.template markDerivativeCoupling<LineLumpedExternalVariables::V1A>();
        signals_.template markDerivativeCoupling<LineLumpedExternalVariables::V1B>();
        signals_.template markDerivativeCoupling<LineLumpedExternalVariables::V1C>();
        signals_.template markDerivativeCoupling<LineLumpedExternalVariables::V2A>();
        signals_.template markDerivativeCoupling<LineLumpedExternalVariables::V2B>();
        signals_.template markDerivativeCoupling<LineLumpedExternalVariables::V2C>();
      }

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
        error_count += z_->verify();
        error_count += y1_->verify();
        error_count += y2_->verify();

        bool matrices_nonzero = hasShuntCapacitance();
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            if (Rp_[n][k] != 0.0 || Lp_[n][k] != 0.0 || Gp_[n][k] != 0.0)
            {
              matrices_nonzero = true;
            }
          }
        }
        if (matrices_nonzero)
        {
          Log::error() << "LineLumped: the rational submodels exclude the "
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
    int LineLumped<scalar_type, index_type>::initialize()
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      for (IdxT j = 0; j < this->ownSize(); ++j)
      {
        y[j]  = 0.0;
        yp[j] = 0.0;
      }

      if (z_.has_value())
      {
        this->initializeSubmodels();
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = true;
      tag_[1] = true;
      tag_[2] = true;
      tag_[3] = false;
      tag_[4] = false;
      tag_[5] = false;
      tag_[6] = false;
      tag_[7] = false;
      tag_[8] = false;

      if (z_.has_value())
      {
        this->tagDifferentiableSubmodels();
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
    int LineLumped<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LineLumped<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y,
        const ScalarT* yp,
        const ScalarT* y_ext,
        const ScalarT* yp_ext,
        ScalarT*       f)
    {
      /* Read variables */
      const ScalarT i12a  = y[0];
      const ScalarT i12b  = y[1];
      const ScalarT i12c  = y[2];
      const ScalarT ish1a = y[3];
      const ScalarT ish1b = y[4];
      const ScalarT ish1c = y[5];
      const ScalarT ish2a = y[6];
      const ScalarT ish2b = y[7];
      const ScalarT ish2c = y[8];

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

      const ScalarT v1a_dot = yp_ext[0];
      const ScalarT v1b_dot = yp_ext[1];
      const ScalarT v1c_dot = yp_ext[2];
      const ScalarT v2a_dot = yp_ext[3];
      const ScalarT v2b_dot = yp_ext[4];
      const ScalarT v2c_dot = yp_ext[5];

      /* 3 series branch equations; the rational series terms accumulate
         through the submodel when the matrix mask is off */
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

      /* 3 terminal 1 shunt algebraic equations */
      f[3] = rl_on_
                 * (G_[0][0] * v1a + G_[0][1] * v1b + G_[0][2] * v1c
                    + C_[0][0] * v1a_dot + C_[0][1] * v1b_dot + C_[0][2] * v1c_dot)
             + TWO<RealT> * ish1a;
      f[4] = rl_on_
                 * (G_[1][0] * v1a + G_[1][1] * v1b + G_[1][2] * v1c
                    + C_[1][0] * v1a_dot + C_[1][1] * v1b_dot + C_[1][2] * v1c_dot)
             + TWO<RealT> * ish1b;
      f[5] = rl_on_
                 * (G_[2][0] * v1a + G_[2][1] * v1b + G_[2][2] * v1c
                    + C_[2][0] * v1a_dot + C_[2][1] * v1b_dot + C_[2][2] * v1c_dot)
             + TWO<RealT> * ish1c;

      /* 3 terminal 2 shunt algebraic equations */
      f[6] = rl_on_
                 * (G_[0][0] * v2a + G_[0][1] * v2b + G_[0][2] * v2c
                    + C_[0][0] * v2a_dot + C_[0][1] * v2b_dot + C_[0][2] * v2c_dot)
             + TWO<RealT> * ish2a;
      f[7] = rl_on_
                 * (G_[1][0] * v2a + G_[1][1] * v2b + G_[1][2] * v2c
                    + C_[1][0] * v2a_dot + C_[1][1] * v2b_dot + C_[1][2] * v2c_dot)
             + TWO<RealT> * ish2b;
      f[8] = rl_on_
                 * (G_[2][0] * v2a + G_[2][1] * v2b + G_[2][2] * v2c
                    + C_[2][0] * v2a_dot + C_[2][1] * v2b_dot + C_[2][2] * v2c_dot)
             + TWO<RealT> * ish2c;

      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LineLumped<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT i12a  = y[0];
      const ScalarT i12b  = y[1];
      const ScalarT i12c  = y[2];
      const ScalarT ish1a = y[3];
      const ScalarT ish1b = y[4];
      const ScalarT ish1c = y[5];
      const ScalarT ish2a = y[6];
      const ScalarT ish2b = y[7];
      const ScalarT ish2c = y[8];

      f_ext[0] = ish1a - i12a;
      f_ext[1] = ish1b - i12b;
      f_ext[2] = ish1c - i12c;
      f_ext[3] = ish2a + i12a;
      f_ext[4] = ish2b + i12b;
      f_ext[5] = ish2c + i12c;

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
      this->evaluateSubmodelInternalResiduals();
      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief External residual contributions to the terminal buses.
     *
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::evaluateExternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f_ext_.data());
      this->scatterExternalResidual();
      this->evaluateSubmodelExternalResiduals();

      return 0;
    }

    /**
     * @brief Residual contribution of the line is pushed to the buses.
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
      monitor_->set(Variable::i_sh1a, [this]
                    { return y_.getData()[3]; });
      monitor_->set(Variable::i_sh1b, [this]
                    { return y_.getData()[4]; });
      monitor_->set(Variable::i_sh1c, [this]
                    { return y_.getData()[5]; });
      monitor_->set(Variable::i_sh2a, [this]
                    { return y_.getData()[6]; });
      monitor_->set(Variable::i_sh2b, [this]
                    { return y_.getData()[7]; });
      monitor_->set(Variable::i_sh2c, [this]
                    { return y_.getData()[8]; });
    }

  } // namespace EMT
} // namespace GridKit

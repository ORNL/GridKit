#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>

#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSource.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase sinusoidal voltage source
     *
     * System sizes:
     * - Number of equations = 6
     * - Number of independent variables = 6
     */
    template <typename scalar_type, typename index_type>
    VoltageSource<scalar_type, index_type>::VoltageSource()
    {
      size_ = 6;
    }

    template <typename scalar_type, typename index_type>
    VoltageSource<scalar_type, index_type>::VoltageSource(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      size_ = 6;
      if (data.Y.has_value())
      {
        own_size_ = 6;
        yfit_.emplace(*data.Y, ONE<RealT>);
        this->registerSubmodel(&*yfit_);
        size_   = own_size_ + yfit_->size();
        rl_on_  = ZERO<RealT>;
        fit_on_ = ONE<RealT>;

        fit_ey_nonzero_ = false;
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            if (data.Y->E[n][k] != 0.0)
            {
              fit_ey_nonzero_ = true;
            }
          }
        }
      }
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    VoltageSource<scalar_type, index_type>::~VoltageSource()
    {
    }

    /**
     * @brief Read model parameters from the data object
     */
    template <typename scalar_type, typename index_type>
    void VoltageSource<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::E))
      {
        E_ = std::get<ABCVector<RealT>>(data.parameters.at(Parameter::E));
      }

      if (data.parameters.contains(Parameter::phi))
      {
        phi_ = std::get<ABCVector<RealT>>(data.parameters.at(Parameter::phi));
      }

      if (data.parameters.contains(Parameter::omega))
      {
        omega_ = std::get<RealT>(data.parameters.at(Parameter::omega));
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
    int VoltageSource<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local storage and registers coupling signals.
     */
    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::allocate()
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
      this->bindPort(u_port_, 3);
      if (yfit_.has_value())
      {
        yfit_->attachInput(&u_port_);
        yfit_->attachOutput(signals_.template getAttachedSignal<VoltageSourceExternalVariables::VA>(),
                            signals_.template getAttachedSignal<VoltageSourceExternalVariables::VB>(),
                            signals_.template getAttachedSignal<VoltageSourceExternalVariables::VC>());
        this->allocateSubmodels();
      }

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize coupling data
      this->allocateExternalVectors(static_cast<IdxT>(VoltageSourceExternalVariables::MAXIMUM), 3);
      signals_.registerExternalVariableSignals(*this);
      this->setExternalResidualSignal(0, signals_.template getAttachedSignal<VoltageSourceExternalVariables::VA>());
      this->setExternalResidualSignal(1, signals_.template getAttachedSignal<VoltageSourceExternalVariables::VB>());
      this->setExternalResidualSignal(2, signals_.template getAttachedSignal<VoltageSourceExternalVariables::VC>());

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     *
     * @return Number of model configuration errors found
     */
    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::verify() const
    {
      int error_count = 0;

      if (!signals_.template isAttached<VoltageSourceExternalVariables::VA>()
          || !signals_.template isAttached<VoltageSourceExternalVariables::VB>()
          || !signals_.template isAttached<VoltageSourceExternalVariables::VC>())
      {
        Log::error() << "VoltageSource: the bus voltage port is not attached\n";
        ++error_count;
      }

      if (omega_ <= 0.0)
      {
        Log::error() << "VoltageSource: the source angular frequency must be positive\n";
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
          Log::error() << "VoltageSource: the rational admittance excludes "
                          "the series matrices\n";
          ++error_count;
        }

        if (fit_ey_nonzero_)
        {
          Log::error() << "VoltageSource: the rational admittance linear "
                          "coefficient must be zero, because the branch "
                          "voltage is algebraic\n";
          ++error_count;
        }
      }

      return error_count;
    }

    /**
     * Initialization of the voltage source model
     *
     */
    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::initialize()
    {
      const RealT sqrt2 = std::numbers::sqrt2_v<RealT>;

      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0] = sqrt2 * E_[0] * std::cos(phi_[0]);
      y[1] = sqrt2 * E_[1] * std::cos(phi_[1]);
      y[2] = sqrt2 * E_[2] * std::cos(phi_[2]);
      y[3] = 0.0;
      y[4] = 0.0;
      y[5] = 0.0;

      yp[0] = 0.0;
      yp[1] = 0.0;
      yp[2] = 0.0;
      yp[3] = 0.0;
      yp[4] = 0.0;
      yp[5] = 0.0;

      if (yfit_.has_value())
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
    int VoltageSource<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = false;
      tag_[1] = false;
      tag_[2] = false;

      // The branch rows carry the differential series current, or the
      // algebraic branch voltage when the rational admittance is present
      bool differential = true;
      if (yfit_.has_value())
      {
        differential = false;
      }
      tag_[3] = differential;
      tag_[4] = differential;
      tag_[5] = differential;

      if (yfit_.has_value())
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
    int VoltageSource<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int VoltageSource<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      /* Read variables */
      const ScalarT ea = y[0];
      const ScalarT eb = y[1];
      const ScalarT ec = y[2];
      const ScalarT ia = y[3];
      const ScalarT ib = y[4];
      const ScalarT ic = y[5];

      /* Read derivatives */
      const ScalarT ia_dot = yp[3];
      const ScalarT ib_dot = yp[4];
      const ScalarT ic_dot = yp[5];

      // Set coupling variable aliases
      const ScalarT va = y_ext[0];
      const ScalarT vb = y_ext[1];
      const ScalarT vc = y_ext[2];

      const RealT sqrt2 = std::numbers::sqrt2_v<RealT>;

      /* 3 source voltage algebraic equations */
      f[0] = ea - sqrt2 * E_[0] * std::cos(omega_ * time_ + phi_[0]);
      f[1] = eb - sqrt2 * E_[1] * std::cos(omega_ * time_ + phi_[1]);
      f[2] = ec - sqrt2 * E_[2] * std::cos(omega_ * time_ + phi_[2]);

      /* 3 series branch equations, series-matrix or rational-admittance
         form; the fit form defines the algebraic branch voltage u read by
         the rational admittance */
      f[3] = rl_on_
                 * (Rs_[0][0] * ia + Rs_[0][1] * ib + Rs_[0][2] * ic
                    + Ls_[0][0] * ia_dot + Ls_[0][1] * ib_dot + Ls_[0][2] * ic_dot)
             + fit_on_ * ia + va - ea;
      f[4] = rl_on_
                 * (Rs_[1][0] * ia + Rs_[1][1] * ib + Rs_[1][2] * ic
                    + Ls_[1][0] * ia_dot + Ls_[1][1] * ib_dot + Ls_[1][2] * ic_dot)
             + fit_on_ * ib + vb - eb;
      f[5] = rl_on_
                 * (Rs_[2][0] * ia + Rs_[2][1] * ib + Rs_[2][2] * ic
                    + Ls_[2][0] * ia_dot + Ls_[2][1] * ib_dot + Ls_[2][2] * ic_dot)
             + fit_on_ * ic + vc - ec;

      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int VoltageSource<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT ia = y[3];
      const ScalarT ib = y[4];
      const ScalarT ic = y[5];

      f_ext[0] = rl_on_ * ia;
      f_ext[1] = rl_on_ * ib;
      f_ext[2] = rl_on_ * ic;

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::evaluateInternalResidual()
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
     * @brief External residual contributions to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::evaluateExternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f_ext_.data());
      this->scatterExternalResidual();
      this->evaluateSubmodelExternalResiduals();

      return 0;
    }

    /**
     * @brief Residual contribution of the source is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* VoltageSource<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void VoltageSource<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::ea, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::eb, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::ec, [this]
                    { return y_.getData()[2]; });
      monitor_->set(Variable::ia, [this]
                    { return y_.getData()[3]; });
      monitor_->set(Variable::ib, [this]
                    { return y_.getData()[4]; });
      monitor_->set(Variable::ic, [this]
                    { return y_.getData()[5]; });
    }

  } // namespace EMT
} // namespace GridKit

#pragma once

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Load/Load.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a pi-model load
     *
     * System sizes:
     * - Number of equations = 2
     * - Number of independent variables = 2
     */
    template <typename scalar_type, typename index_type>
    Load<scalar_type, index_type>::Load(BusT* bus)
      : bus_(bus)
    {
      size_ = 2;
    }

    template <typename scalar_type, typename index_type>
    Load<scalar_type, index_type>::Load(BusT* bus,
                                        RealT R,
                                        RealT X)
      : bus_(bus),
        R_(R),
        X_(X)
    {
      size_ = 2;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    Load<scalar_type, index_type>::Load(BusT*             bus,
                                        const ModelDataT& data)
      : bus_(bus)
    {
      if (data.parameters.contains(ModelDataT::Parameters::R))
      {
        R_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::R));
      }

      if (data.parameters.contains(ModelDataT::Parameters::X))
      {
        X_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::X));
      }

      // using Variable = typename ModelDataT::MonitorableVariables;
      // monitor_->set(Variable::p, [this] { return ?; });
      // monitor_->set(Variable::q, [this] { return ?; });

      size_ = 2;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    Load<scalar_type, index_type>::~Load()
    {
      // std::cout << "Destroy Load..." << std::endl;
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::allocate()
    {
      // std::cout << "Allocate Load..." << std::endl;

      auto size = static_cast<size_t>(size_); // avoid compiler warnings
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Resize coupling data
      wb_.resize(2);
      h_.resize(2);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      return 0;
    }

    /**
     * Initialization of the load model
     *
     */
    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::initialize()
    {
      ScalarT vr = Vr();
      ScalarT vi = Vi();
      ScalarT ir = -(g_ * vr - b_ * vi);
      ScalarT ii = -(b_ * vr + g_ * vi);

      y_[0] = ir;
      y_[1] = ii;

      yp_[0] = 0.0;
      yp_[1] = 0.0;

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = false;
      tag_[1] = false;

      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Load<scalar_type, index_type>::evaluateBusResidual(
        ScalarT*                  y,
        [[maybe_unused]] ScalarT* yp,
        [[maybe_unused]] ScalarT* wb,
        ScalarT*                  h)
    {
      ScalarT Ir = y[0];
      ScalarT Ii = y[1];
      h[0]       = Ir;
      h[1]       = Ii;

      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Load<scalar_type, index_type>::evaluateInternalResidual(
        ScalarT*                  y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  f)
    {
      ScalarT Vr = wb[0];
      ScalarT Vi = wb[1];
      ScalarT Ir = y[0];
      ScalarT Ii = y[1];
      f[0]       = Ir + g_ * Vr - b_ * Vi;
      f[1]       = Ii + b_ * Vr + g_ * Vi;

      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::evaluateResidual()
    {
      wb_[0] = Vr();
      wb_[1] = Vi();
      evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), f_.data());
      evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());
      Ir() += h_[0];
      Ii() += h_[1];

      return 0;
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename scalar_type, typename index_type>
    void Load<scalar_type, index_type>::setDerivedParams()
    {
      b_ = -X_ / (R_ * R_ + X_ * X_);
      g_ = R_ / (R_ * R_ + X_ * X_);
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Load<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

  } // namespace PhasorDynamics
} // namespace GridKit

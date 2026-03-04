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
     * - Number of equations = 0
     * - Number of independent variables = 0
     */
    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type* bus)
      : bus_(bus)
    {
      size_ = 0;
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type* bus,
                              RealT     R,
                              RealT     X)
      : bus_(bus),
        R_(R),
        X_(X)
    {
      size_ = 0;
      setDerivedParams();
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type*              bus,
                              const model_data_type& data)
      : bus_(bus)
    {
      if (data.parameters.contains(model_data_type::Parameters::R))
      {
        R_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::R));
      }

      if (data.parameters.contains(model_data_type::Parameters::X))
      {
        X_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::X));
      }

      // using Variable = typename model_data_type::MonitorableVariables;
      // monitor_->set(Variable::p, [this] { return ?; });
      // monitor_->set(Variable::q, [this] { return ?; });

      size_ = 0;
      setDerivedParams();
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::~Load()
    {
      // std::cout << "Destroy Load..." << std::endl;
    }

    /**
     * @brief Set the component ID
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::allocate()
    {
      // std::cout << "Allocate Load..." << std::endl;

      wb_.resize(2);
      h_.resize(2);

      return 0;
    }

    /**
     * Initialization of the load model
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::tagDifferentiable()
    {
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
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int Load<ScalarT, IdxT>::evaluateBusResidual(
        [[maybe_unused]] ScalarT* y, [[maybe_unused]] ScalarT* yp, ScalarT* wb, ScalarT* h)
    {
      ScalarT Vr = wb[0];
      ScalarT Vi = wb[1];
      ScalarT Ir = -g_ * Vr + b_ * Vi;
      ScalarT Ii = -b_ * Vr - g_ * Vi;
      h[0]       = Ir;
      h[1]       = Ii;

      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateResidual()
    {
      wb_[0] = Vr();
      wb_[1] = Vi();
      evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());
      Ir() += h_[0];
      Ii() += h_[1];

      return 0;
    }

    /**
     * @brief Derived parameters
     *
     */
    template <class ScalarT, typename IdxT>
    void Load<ScalarT, IdxT>::setDerivedParams()
    {
      b_ = -X_ / (R_ * R_ + X_ * X_);
      g_ = R_ / (R_ * R_ + X_ * X_);
    }

    template <class ScalarT, typename IdxT>
    const Model::VariableMonitorBase* Load<ScalarT, IdxT>::getMonitor() const
    {
      return monitor_.get();
    }

  } // namespace PhasorDynamics
} // namespace GridKit

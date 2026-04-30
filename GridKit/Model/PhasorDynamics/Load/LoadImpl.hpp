#pragma once

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElementImpl.hpp>
#include <GridKit/Model/PhasorDynamics/Load/Load.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadData.hpp>

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
    template <typename ScalarP, typename IdxP>
    Load<ScalarP, IdxP>::Load(BusT* bus)
      : bus_(bus)
    {
      size_ = 0;
    }

    template <typename ScalarP, typename IdxP>
    Load<ScalarP, IdxP>::Load(BusT* bus, RealT R, RealT X)
      : bus_(bus),
        R_(R),
        X_(X)
    {
      size_ = 0;
      setDerivedParams();
    }

    template <typename ScalarP, typename IdxP>
    Load<ScalarP, IdxP>::Load(BusT* bus, const ModelDataT& data)
      : ConnectedElement<Load>(data),
        bus_(bus)
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

      size_ = 0;
      setDerivedParams();
    }

    template <typename ScalarP, typename IdxP>
    Load<ScalarP, IdxP>::~Load()
    {
      // std::cout << "Destroy Load..." << std::endl;
    }

    /**
     * @brief Set the component ID
     */
    template <typename ScalarP, typename IdxP>
    int Load<ScalarP, IdxP>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename ScalarP, typename IdxP>
    int Load<ScalarP, IdxP>::allocate()
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
    template <typename ScalarP, typename IdxP>
    int Load<ScalarP, IdxP>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename ScalarP, typename IdxP>
    int Load<ScalarP, IdxP>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <typename ScalarP, typename IdxP>
    __attribute__((always_inline)) int Load<ScalarP, IdxP>::evaluateBusResidual(
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
    template <typename ScalarP, typename IdxP>
    int Load<ScalarP, IdxP>::evaluateResidual()
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
    template <typename ScalarP, typename IdxP>
    void Load<ScalarP, IdxP>::setDerivedParams()
    {
      b_ = -X_ / (R_ * R_ + X_ * X_);
      g_ = R_ / (R_ * R_ + X_ * X_);
    }

  } // namespace PhasorDynamics
} // namespace GridKit

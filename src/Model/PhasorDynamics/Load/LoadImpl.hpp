#pragma once

#include <cmath>
#include <iostream>

#include "Load.hpp"
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Load/LoadData.hpp>

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
                              real_type R,
                              real_type X)
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
        R_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::R));
      }

      if (data.parameters.contains(model_data_type::Parameters::X))
      {
        X_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::X));
      }

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

      w_.resize(2);
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
     * @brief Bus residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int Load<ScalarT, IdxT>::evaluateBusResidual(
        [[maybe_unused]] ScalarT* y, [[maybe_unused]] ScalarT* yp, ScalarT* w, ScalarT* h)
    {
      ScalarT Vr = w[0];
      ScalarT Vi = w[1];
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
      w_[0] = Vr();
      w_[1] = Vi();
      evaluateBusResidual(y_.data(), yp_.data(), w_.data(), h_.data());
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

  } // namespace PhasorDynamics
} // namespace GridKit

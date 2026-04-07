#pragma once

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/LoadZIP/LoadZIP.hpp>
#include <GridKit/Model/PhasorDynamics/LoadZIP/LoadZIPData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a zip load
     *
     * System sizes:
     * - Number of equations = 2
     * - Number of independent variables = 2
     */
    template <class ScalarT, typename IdxT>
    LoadZIP<ScalarT, IdxT>::LoadZIP(bus_type* bus)
      : bus_(bus)
    {
      size_ = 2;
    }

    template <class ScalarT, typename IdxT>
    LoadZIP<ScalarT, IdxT>::LoadZIP(bus_type* bus, RealT P0, RealT Q0, RealT V0, RealT alphaI, RealT alphaP)
      : bus_(bus),
        P0_(P0),
        Q0_(Q0),
        V0_(V0),
        alphaI_(alphaI),
        alphaP_(alphaP)
    {
      size_ = 2;
      setDerivedParams();
    }

    template <class ScalarT, typename IdxT>
    LoadZIP<ScalarT, IdxT>::LoadZIP(bus_type*              bus,
                                    const model_data_type& data)
      : bus_(bus)
    {
      if (data.parameters.contains(model_data_type::Parameters::P0))
      {
        P0_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::P0));
      }

      if (data.parameters.contains(model_data_type::Parameters::Q0))
      {
        Q0_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Q0));
      }

      if (data.parameters.contains(model_data_type::Parameters::V0))
      {
        V0_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::V0));
      }

      if (data.parameters.contains(model_data_type::Parameters::alphaI))
      {
        alphaI_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::alphaI));
      }

      if (data.parameters.contains(model_data_type::Parameters::alphaP))
      {
        alphaP_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::alphaP));
      }

      // using Variable = typename model_data_type::MonitorableVariables;
      // monitor_->set(Variable::p, [this] { return ?; });
      // monitor_->set(Variable::q, [this] { return ?; });

      size_ = 2;
      setDerivedParams();
    }

    template <class ScalarT, typename IdxT>
    LoadZIP<ScalarT, IdxT>::~LoadZIP()
    {
      // std::cout << "Destroy LoadZIP..." << std::endl;
    }

    /**
     * @brief Set the component ID
     */
    template <class ScalarT, typename IdxT>
    int LoadZIP<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int LoadZIP<ScalarT, IdxT>::allocate()
    {
      // std::cout << "Allocate Load..." << std::endl;

      auto size = static_cast<size_t>(size_); // avoid compiler warnings
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      wb_.resize(2);

      ws_.resize(1);
      ws_indices_.resize(1);
      ws_[0]         = 0.0;
      ws_indices_[0] = INVALID_INDEX<IdxT>;

      h_.resize(2);

      return 0;
    }

    /**
     * Initialization of the load model
     *
     */
    template <class ScalarT, typename IdxT>
    int LoadZIP<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int LoadZIP<ScalarT, IdxT>::tagDifferentiable()
    {
      tag_[0] = false;  
      tag_[1] = false;
      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int LoadZIP<ScalarT, IdxT>::evaluateBusResidual(
        [[maybe_unused]] ScalarT* y, [[maybe_unused]] ScalarT* yp, [[maybe_unused]] ScalarT* wb, ScalarT* h)
    {
      ScalarT Ir = y[0];
      ScalarT Ii = y[1];
      h[0]       = Ir;
      h[1]       = Ii;

      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <class ScalarT, typename IdxT>
    int LoadZIP<ScalarT, IdxT>::evaluateResidual()
    {
      wb_[0] = Vr();
      wb_[1] = Vi();
      evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
      evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());
      Ir() += h_[0];
      Ii() += h_[1];

      return 0;
    }

    /**
     * @brief Internal residual
     * 
     */
    template <class ScalarT, typename IdxT>
    int LoadZIP<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
      [[maybe_unused]] ScalarT* yp, ScalarT* wb, [[maybe_unused]] ScalarT* ws, ScalarT* f)
    {
      ScalarT one{1.0};
      ScalarT Vr = wb[0];
      ScalarT Vi = wb[1];
      ScalarT Ir = y[0];
      ScalarT Ii = y[1];
      ScalarT Vm2 = Vr*Vr + Vi*Vi;
      ScalarT Vm = std::sqrt(Vm2);
      ScalarT ifrac = (one / (V0_*V0_) * (one - alphaI_ - alphaP_)
        + one / (V0_ * Vm) * alphaI_ + one / Vm2 * alphaP_);
      f[0] = Ir + (P0_ * Vr + Q0_ * Vi) * ifrac;
      f[1] = Ii + (P0_ * Vi - Q0_ * Vr) * ifrac;
      return 0;
    }

    /**
     * @brief Derived parameters
     *
     */
    template <class ScalarT, typename IdxT>
    void LoadZIP<ScalarT, IdxT>::setDerivedParams()
    {
      return;
    }

    template <class ScalarT, typename IdxT>
    const Model::VariableMonitorBase* LoadZIP<ScalarT, IdxT>::getMonitor() const
    {
      return monitor_.get();
    }

  } // namespace PhasorDynamics
} // namespace GridKit

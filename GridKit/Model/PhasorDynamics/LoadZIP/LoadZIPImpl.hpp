#pragma once

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElementImpl.hpp>
#include <GridKit/Model/PhasorDynamics/LoadZIP/LoadZIP.hpp>
#include <GridKit/Model/PhasorDynamics/LoadZIP/LoadZIPData.hpp>

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
    template <typename ScalarP, typename IdxP>
    LoadZIP<ScalarP, IdxP>::LoadZIP(BusT* bus)
      : bus_(bus)
    {
      size_ = 2;
    }

    template <typename ScalarP, typename IdxP>
    LoadZIP<ScalarP, IdxP>::LoadZIP(BusT* bus, RealT P0, RealT Q0, RealT V0, RealT alphaI, RealT alphaP)
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

    template <typename ScalarP, typename IdxP>
    LoadZIP<ScalarP, IdxP>::LoadZIP(BusT* bus, const ModelDataT& data)
      : ConnectedElement<LoadZIP>(data),
        bus_(bus)
    {
      if (data.parameters.contains(ModelDataT::Parameters::P0))
      {
        P0_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::P0));
      }

      if (data.parameters.contains(ModelDataT::Parameters::Q0))
      {
        Q0_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Q0));
      }

      if (data.parameters.contains(ModelDataT::Parameters::V0))
      {
        V0_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::V0));
      }

      if (data.parameters.contains(ModelDataT::Parameters::alphaI))
      {
        alphaI_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::alphaI));
      }

      if (data.parameters.contains(ModelDataT::Parameters::alphaP))
      {
        alphaP_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::alphaP));
      }

      // using Variable = typename ModelDataT::MonitorableVariables;
      // monitor_->set(Variable::p, [this] { return ?; });
      // monitor_->set(Variable::q, [this] { return ?; });

      size_ = 2;
      setDerivedParams();
    }

    template <typename ScalarP, typename IdxP>
    LoadZIP<ScalarP, IdxP>::~LoadZIP()
    {
      // std::cout << "Destroy LoadZIP..." << std::endl;
    }

    /**
     * @brief Set the component ID
     */
    template <typename ScalarP, typename IdxP>
    int LoadZIP<ScalarP, IdxP>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename ScalarP, typename IdxP>
    int LoadZIP<ScalarP, IdxP>::allocate()
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
    template <typename ScalarP, typename IdxP>
    int LoadZIP<ScalarP, IdxP>::initialize()
    {
      ScalarT vr    = Vr();
      ScalarT vi    = Vi();
      ScalarT Vm2   = vr * vr + vi * vi;
      ScalarT Vm    = std::sqrt(Vm2);
      ScalarT ifrac = (ONE<RealT> - alphaI_ - alphaP_) / (V0_ * V0_)
                      + alphaI_ / (V0_ * Vm)
                      + alphaP_ / Vm2;
      ScalarT ir = -(P0_ * vr + Q0_ * vi) * ifrac;
      ScalarT ii = -(P0_ * vi - Q0_ * vr) * ifrac;
      y_[0]      = ir;
      y_[1]      = ii;

      yp_[0] = 0.0;
      yp_[1] = 0.0;
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename ScalarP, typename IdxP>
    int LoadZIP<ScalarP, IdxP>::tagDifferentiable()
    {
      tag_[0] = false;
      tag_[1] = false;
      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <typename ScalarP, typename IdxP>
    __attribute__((always_inline)) int LoadZIP<ScalarP, IdxP>::evaluateBusResidual(
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
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <typename ScalarP, typename IdxP>
    int LoadZIP<ScalarP, IdxP>::evaluateResidual()
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
     * @brief Internal residual
     *
     */
    template <typename ScalarP, typename IdxP>
    __attribute__((always_inline)) int LoadZIP<ScalarP, IdxP>::evaluateInternalResidual(
        ScalarT*                  y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  f)
    {
      ScalarT Vr    = wb[0];
      ScalarT Vi    = wb[1];
      ScalarT Ir    = y[0];
      ScalarT Ii    = y[1];
      ScalarT Vm2   = Vr * Vr + Vi * Vi;
      ScalarT Vm    = std::sqrt(Vm2);
      ScalarT ifrac = (ONE<RealT> - alphaI_ - alphaP_) / (V0_ * V0_)
                      + alphaI_ / (V0_ * Vm)
                      + alphaP_ / Vm2;
      f[0] = Ir + (P0_ * Vr + Q0_ * Vi) * ifrac;
      f[1] = Ii + (P0_ * Vi - Q0_ * Vr) * ifrac;
      return 0;
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename ScalarP, typename IdxP>
    void LoadZIP<ScalarP, IdxP>::setDerivedParams()
    {
      return;
    }

  } // namespace PhasorDynamics
} // namespace GridKit

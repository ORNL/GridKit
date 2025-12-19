#pragma once

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Constructor for a bus fault
     *
     * Model sizes:
     * - Number of equations = 0
     * - Number of independent variables = 0
     */
    template <class ScalarT, typename IdxT>
    BusFault<ScalarT, IdxT>::BusFault(bus_type* bus)
      : bus_(bus), R_(0), X_(0.01), status_(0), bus_id_(0)
    {
      (void) bus_id_;
      size_ = 0;
      setDerivedParams();
    }

    /**
     * @brief Construct a new BusFault
     *
     * @tparam ScalarT - scalar type
     * @tparam IdxT    - matrix/vector index type
     * @param bus1 - pointer to bus-1
     * @param bus2 - pointer to bus-2
     * @param R - line series resistance
     * @param X - line series reactance
     * @param G - line shunt conductance
     * @param B - line shunt charging
     */
    template <class ScalarT, typename IdxT>
    BusFault<ScalarT, IdxT>::BusFault(bus_type* bus, RealT R, RealT X, int status)
      : bus_(bus), R_(R), X_(X), status_(status), bus_id_(0)
    {
      size_ = 0;
      setDerivedParams();
    }

    /**
     * @brief Construct a new BusFault
     *
     * @tparam ScalarT - scalar type
     * @tparam IdxT    - matrix/vector index type
     * @param bus1 - pointer to bus-1
     * @param bus2 - pointer to bus-2
     */
    template <class ScalarT, typename IdxT>
    BusFault<ScalarT, IdxT>::BusFault(bus_type* bus, const DataT& data)
      : bus_(bus)
    {
      if (data.parameters.contains(DataT::Parameters::R))
      {
        R_ = std::get<RealT>(data.parameters.at(DataT::Parameters::R));
      }

      if (data.parameters.contains(DataT::Parameters::X))
      {
        X_ = std::get<RealT>(data.parameters.at(DataT::Parameters::X));
      }

      if (data.parameters.contains(DataT::Parameters::state0))
      {
        status_ = std::get<bool>(data.parameters.at(DataT::Parameters::state0));
      }

      if (data.ports.contains(DataT::Ports::bus))
      {
        bus_id_ = data.ports.at(DataT::Ports::bus);
      }

      using Variable = typename DataT::MonitorableVariables;
      // monitor_.set(Variable::state, [this]
      //              { return status_; });
      monitor_.set(Variable::ir, [this]
                   { return Ir(); });
      monitor_.set(Variable::ii, [this]
                   { return Ii(); });

      size_ = 0;
      setDerivedParams();
    }

    /**
     * @brief Set the component ID
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::allocate()
    {
      // std::cout << "Allocate BusFault..." << std::endl;

      wb_.resize(2);
      h_.resize(2);

      return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int BusFault<ScalarT, IdxT>::evaluateBusResidual(
        [[maybe_unused]] ScalarT* y, [[maybe_unused]] ScalarT* yp, ScalarT* wb, ScalarT* h)
    {
      ScalarT Vr = wb[0];
      ScalarT Vi = wb[1];
      ScalarT Ir = -Vr * G_ + Vi * B_;
      ScalarT Ii = -Vr * B_ - Vi * G_;
      h[0]       = Ir;
      h[1]       = Ii;

      return 0;
    }

    /**
     * \brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     *
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::evaluateResidual()
    {
      if (status_)
      {
        wb_[0] = Vr();
        wb_[1] = Vi();
        evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());
        Ir() += h_[0];
        Ii() += h_[1];
      }
      return 0;
    }

    /**
     * @brief Derived parameters
     *
     */
    template <class ScalarT, typename IdxT>
    void BusFault<ScalarT, IdxT>::setDerivedParams()
    {
      B_ = -X_ / (X_ * X_ + R_ * R_);
      G_ = R_ / (X_ * X_ + R_ * R_);
    }
  } // namespace PhasorDynamics
} // namespace GridKit

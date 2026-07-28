#pragma once

#include <algorithm>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a constant-impedance load
     *
     * Model size:
     * - Number of equations = 0
     * - Number of internal variables = 0
     */
    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(BusT* bus)
      : bus_(bus)
    {
      size_ = 0;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(BusT* bus,
                                          RealT R,
                                          RealT X)
      : bus_(bus),
        R_(R),
        X_(X)
    {
      size_ = 0;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(BusT*             bus,
                                          const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::R))
      {
        R_ = std::get<RealT>(data.parameters.at(Parameter::R));
      }

      if (data.parameters.contains(Parameter::X))
      {
        X_ = std::get<RealT>(data.parameters.at(Parameter::X));
      }

      size_ = 0;
      setDerivedParams();
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::~LoadZ()
    {
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
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      allocated_ = true;
      return 0;
    }

    /**
     * Initialization of the load model
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::tagDifferentiable()
    {
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
    int LoadZ<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateResidual()
    {
      const ScalarT Vr = this->Vr();
      const ScalarT Vi = this->Vi();

      Ir() -= g_ * Vr - b_ * Vi;
      Ii() -= b_ * Vr + g_ * Vi;

      return 0;
    }

    /**
     * @brief A constant impedance load is one constant admittance stamp.
     *
     * The sign is negative because the load draws current from the bus, matching
     * the subtraction in evaluateResidual().
     *
     */
    template <typename scalar_type, typename index_type>
    index_type LoadZ<scalar_type, index_type>::admittanceStamps(
        typename Component<ScalarT, IdxT>::StampT* out)
    {
      // See Branch::admittanceStamps for why an infinite bus terminal opts out.
      if (bus_->size() == 0)
      {
        return 0;
      }

      if (out != nullptr)
      {
        out[0] = {bus_->getResidualIndices()[0], bus_->getVariableIndices()[0], -g_, -b_};
      }

      return 1;
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename scalar_type, typename index_type>
    void LoadZ<scalar_type, index_type>::setDerivedParams()
    {
      b_ = -X_ / (R_ * R_ + X_ * X_);
      g_ = R_ / (R_ * R_ + X_ * X_);
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

      monitor_->set(Variable::p, [this]
                    {
                      const ScalarT Vr = this->Vr();
                      const ScalarT Vi = this->Vi();
                      return -g_ * (Vr * Vr + Vi * Vi); });
      monitor_->set(Variable::q, [this]
                    {
                      const ScalarT Vr = this->Vr();
                      const ScalarT Vi = this->Vi();
                      return b_ * (Vr * Vr + Vi * Vi); });
    }

  } // namespace PhasorDynamics
} // namespace GridKit

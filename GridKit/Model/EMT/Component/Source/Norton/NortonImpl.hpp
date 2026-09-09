/**
 * @file NortonImpl.hpp
 * @brief Implementation of the EMT Norton model.
 */
#pragma once

#include <cmath>
#include <stdexcept>

#include <GridKit/Model/EMT/Component/Source/Norton/Norton.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Norton<scalar_type, index_type>::Norton(const YDataT& Y, PhaseSignals voltage, RealT scale)
      : voltage_(voltage), admittance_(Y, scale)
    {
      if (Y.rows != 3 || Y.cols != 3 || Y.validate() != 0 || !std::isfinite(scale))
        throw std::invalid_argument("Norton: expected a finite three-phase admittance");
      equation_size_ = 3;
      size_          = 3 + admittance_.size();
      this->addOperator(&admittance_);
      for (auto& signal : shunt_)
        signal.claimProducer();
      admittance_.attachInput(voltage_[0], voltage_[1], voltage_[2]);
      admittance_.attachOutput(&shunt_[0], &shunt_[1], &shunt_[2]);
    }

    template <typename scalar_type, typename index_type>
    typename Norton<scalar_type, index_type>::SignalT& Norton<scalar_type, index_type>::outputSignal(size_t phase)
    {
      return shunt_.at(phase);
    }

    template <typename scalar_type, typename index_type>
    const typename Norton<scalar_type, index_type>::SignalT& Norton<scalar_type, index_type>::outputSignal(size_t phase) const
    {
      return shunt_.at(phase);
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
        this->allocateVectors(size_);
      tag_.resize(static_cast<size_t>(size_));
      variable_indices_.resize(static_cast<size_t>(size_));
      residual_indices_.resize(static_cast<size_t>(size_));
      for (IdxT p = 0; p < 3; ++p)
        this->bindSignal(shunt_[static_cast<size_t>(p)], p);
      const int status = this->allocateOperators();
      if (status != 0)
        return status;
      this->assignGlobalIndices(0);
      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::verify() const
    {
      return admittance_.verify();
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
    {
      if (!values.empty())
        throw std::invalid_argument("Norton shunt current is initialized from its admittance");
      return initialize();
    }

    template <typename scalar_type, typename index_type>
    typename Component<scalar_type, index_type>::InitializationPortsT Norton<scalar_type, index_type>::initializationPorts()
    {
      return {{voltage_.begin(), voltage_.end()}, {}, {}};
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::initialize()
    {
      const int status = admittance_.initialize();
      if (status != 0)
        return status;
      return initializeShunt();
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::initializeShunt()
    {
      for (size_t p = 0; p < 3; ++p)
      {
        shunt_[p].init(admittance_.output(static_cast<IdxT>(p)));
        shunt_[p].initDerivative(ScalarT{0});
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::initializeSteadyState(RealT omega)
    {
      const int status = admittance_.initializeSteadyState(omega);
      return status == 0 ? initializeShunt() : status;
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
      return this->setAbsoluteToleranceOperators(tolerance);
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Norton<scalar_type, index_type>::evaluateInternalResidual(const ScalarT* y,
                                                                                                 const ScalarT*,
                                                                                                 const ScalarT*,
                                                                                                 const ScalarT*,
                                                                                                 ScalarT* f)
    {
      for (size_t p = 0; p < 3; ++p)
        f[p] = -y[p];
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::evaluateInternalResidual()
    {
      evaluateInternalResidual(y_.getData(), nullptr, nullptr, nullptr, f_.getData());
      const int status = this->evaluateOperatorInternalResiduals();
      f_.setDataUpdated();
      return status;
    }

    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::evaluateResidual()
    {
      const int status = evaluateInternalResidual();
      return status == 0 ? this->evaluateExternalResidual() : status;
    }
  } // namespace EMT
} // namespace GridKit

/**
 * @file FilterImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include <GridKit/Model/EMT/Component/Filter/Filter.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Filter<scalar_type, index_type>::Filter(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      equation_size_ = size_ = static_cast<IdxT>(FilterInternalVariables::MAXIMUM);
      for (size_t n = 0; n < output_.size(); ++n)
      {
        signals_.assignSignal(static_cast<FilterInternalVariables>(n), &output_[n]);
      }
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Filter<scalar_type, index_type>::~Filter() = default;

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      Rs_             = parameter<ABCMatrix<RealT>>(data, Parameter::Rs, Rs_);
      Ls_             = parameter<ABCMatrix<RealT>>(data, Parameter::Ls);
      C_              = parameter<ABCMatrix<RealT>>(data, Parameter::C);
      Rg_             = parameter<ABCMatrix<RealT>>(data, Parameter::Rg, Rg_);
      Lg_             = parameter<ABCMatrix<RealT>>(data, Parameter::Lg);

      // Principal minors characterize positive (semi)definite symmetric
      // three-phase matrices. Normalize to keep the check independent of units.
      auto passive = [](const ABCMatrix<RealT>& matrix, bool strict)
      {
        RealT scale = ZERO<RealT>;
        for (const auto& row : matrix)
        {
          for (const auto value : row)
          {
            if (!std::isfinite(value))
            {
              return false;
            }
            scale = std::max(scale, std::abs(value));
          }
        }
        if (scale == ZERO<RealT>)
        {
          return !strict;
        }
        ABCMatrix<RealT> a{};
        const RealT      tolerance = 3 * std::numeric_limits<RealT>::epsilon();
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            a[n][k] = matrix[n][k] / scale;
            if (std::abs(a[n][k] - matrix[k][n] / scale) > tolerance)
            {
              return false;
            }
          }
        }
        auto positive = [strict, tolerance](RealT value)
        {
          return strict ? value > tolerance : value >= -tolerance;
        };
        for (size_t n = 0; n < 3; ++n)
        {
          if (!positive(a[n][n]))
          {
            return false;
          }
          for (size_t k = n + 1; k < 3; ++k)
          {
            if (!positive(a[n][n] * a[k][k] - a[n][k] * a[k][n]))
            {
              return false;
            }
          }
        }
        return positive(a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
      };
      if (!passive(Rs_, false) || !passive(Rg_, false))
      {
        throw std::invalid_argument("Filter: Rs and Rg must be finite symmetric positive semidefinite matrices");
      }
      if (!passive(Ls_, true) || !passive(C_, true) || !passive(Lg_, true))
      {
        throw std::invalid_argument("Filter: Ls, C and Lg must be finite symmetric positive definite matrices");
      }
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::attachInput(PhaseSignals voltage, PhaseSignals source)
    {
      if (allocated_)
      {
        throw std::logic_error("Filter inputs cannot change after allocation");
      }
      for (size_t p = 0; p < 3; ++p)
      {
        signals_.attachSignal(static_cast<FilterExternalVariables>(p), voltage[p]);
        signals_.attachSignal(static_cast<FilterExternalVariables>(3 + p), source[p]);
      }
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      if (allocated_)
      {
        throw std::logic_error("Filter outputs cannot change after allocation");
      }
      const auto n = static_cast<size_t>(output);
      if (n >= output_.size() || !signal || alias_[n])
      {
        throw std::invalid_argument("Filter: invalid or duplicate output assignment");
      }
      signal->claimProducer();
      alias_[n] = signal;
    }

    template <typename scalar_type, typename index_type>
    typename Filter<scalar_type, index_type>::SignalT& Filter<scalar_type, index_type>::outputSignal(Outputs output)
    {
      return output_.at(static_cast<size_t>(output));
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::allocate()
    {
      using V = FilterExternalVariables;
      if (signals_.attachedSignals({V::VA, V::VB, V::VC, V::EA, V::EB, V::EC}).size() != 6)
      {
        throw std::invalid_argument("Filter: all voltage inputs are required");
      }
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      tag_.resize(static_cast<size_t>(size_));
      variable_indices_.resize(static_cast<size_t>(size_));
      residual_indices_.resize(static_cast<size_t>(equation_size_));
      this->assignGlobalIndices(0);
      this->allocateExternalVectors(static_cast<IdxT>(V::MAXIMUM), 0);
      signals_.registerExternalVariableSignals(*this);
      signals_.bindInternalVariableSignals(*this);
      for (size_t n = 0; n < alias_.size(); ++n)
      {
        if (alias_[n])
        {
          this->bindSignal(*alias_[n], static_cast<IdxT>(n));
        }
      }
      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::verify() const
    {
      using V           = FilterExternalVariables;
      const auto inputs = signals_.attachedSignals({V::VA, V::VB, V::VC, V::EA, V::EB, V::EC});
      int        errors = inputs.size() != 6;
      for (const auto* input : inputs)
      {
        errors += !input->linked();
      }
      return errors;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      for (size_t n = 0; n < output_.size(); ++n)
      {
        output_[n].init(static_cast<ScalarT>(this->outputValue(outputs, static_cast<Outputs>(n), ZERO<RealT>)));
        output_[n].initDerivative(ScalarT{0});
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
    {
      this->template parseInitialOutputs<Filter>(values);
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
    {
      return this->initializeOutputs(*this, values);
    }

    template <typename scalar_type, typename index_type>
    typename Component<scalar_type, index_type>::InitializationPortsT Filter<scalar_type, index_type>::initializationPorts()
    {
      typename Component<ScalarT, IdxT>::InitializationPortsT ports;
      for (const auto output : magic_enum::enum_values<Outputs>())
      {
        if (output != Outputs::SIZE)
        {
          ports.outputs.emplace(std::string(magic_enum::enum_name(output)), &outputSignal(output));
        }
      }
      return ports;
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::prepareInitialization(typename Component<ScalarT, IdxT>::InitialStateT& initial)
    {
      const auto outputs = this->template parseInitialOutputs<Filter>(initial.outputs(*this));
      for (size_t n = 0; n < output_.size(); ++n)
      {
        const auto value = this->outputValue(outputs, static_cast<Outputs>(n), ZERO<RealT>);
        initial.provide(output_[n], value);
        if (alias_[n])
        {
          initial.provide(*alias_[n], value);
        }
      }
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
      return 0;
    }

    /**
     * @brief Internal residual
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Filter<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      /* Read variables */
      const ScalarT ia  = y[0];
      const ScalarT ib  = y[1];
      const ScalarT ic  = y[2];
      const ScalarT voa = y[3];
      const ScalarT vob = y[4];
      const ScalarT voc = y[5];
      const ScalarT iga = y[6];
      const ScalarT igb = y[7];
      const ScalarT igc = y[8];

      /* Read derivatives */
      const ScalarT ia_dot  = yp[0];
      const ScalarT ib_dot  = yp[1];
      const ScalarT ic_dot  = yp[2];
      const ScalarT voa_dot = yp[3];
      const ScalarT vob_dot = yp[4];
      const ScalarT voc_dot = yp[5];
      const ScalarT iga_dot = yp[6];
      const ScalarT igb_dot = yp[7];
      const ScalarT igc_dot = yp[8];

      // Set coupling variable aliases
      const ScalarT va = y_ext[0];
      const ScalarT vb = y_ext[1];
      const ScalarT vc = y_ext[2];
      const ScalarT ea = y_ext[3];
      const ScalarT eb = y_ext[4];
      const ScalarT ec = y_ext[5];

      /* Converter-side inductor equations */
      f[0] = voa - ea
             + (Rs_[0][0] * ia + Ls_[0][0] * ia_dot)
             + (Rs_[0][1] * ib + Ls_[0][1] * ib_dot)
             + (Rs_[0][2] * ic + Ls_[0][2] * ic_dot);
      f[1] = vob - eb
             + (Rs_[1][0] * ia + Ls_[1][0] * ia_dot)
             + (Rs_[1][1] * ib + Ls_[1][1] * ib_dot)
             + (Rs_[1][2] * ic + Ls_[1][2] * ic_dot);
      f[2] = voc - ec
             + (Rs_[2][0] * ia + Ls_[2][0] * ia_dot)
             + (Rs_[2][1] * ib + Ls_[2][1] * ib_dot)
             + (Rs_[2][2] * ic + Ls_[2][2] * ic_dot);

      /* Capacitor equations */
      f[3] = iga - ia + C_[0][0] * voa_dot + C_[0][1] * vob_dot + C_[0][2] * voc_dot;
      f[4] = igb - ib + C_[1][0] * voa_dot + C_[1][1] * vob_dot + C_[1][2] * voc_dot;
      f[5] = igc - ic + C_[2][0] * voa_dot + C_[2][1] * vob_dot + C_[2][2] * voc_dot;

      /* Grid-side inductor equations */
      f[6] = va - voa
             + (Rg_[0][0] * iga + Lg_[0][0] * iga_dot)
             + (Rg_[0][1] * igb + Lg_[0][1] * igb_dot)
             + (Rg_[0][2] * igc + Lg_[0][2] * igc_dot);
      f[7] = vb - vob
             + (Rg_[1][0] * iga + Lg_[1][0] * iga_dot)
             + (Rg_[1][1] * igb + Lg_[1][1] * igb_dot)
             + (Rg_[1][2] * igc + Lg_[1][2] * igc_dot);
      f[8] = vc - voc
             + (Rg_[2][0] * iga + Lg_[2][0] * iga_dot)
             + (Rg_[2][1] * igb + Lg_[2][1] * igb_dot)
             + (Rg_[2][2] * igc + Lg_[2][2] * igc_dot);

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();
      evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
      f_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    void Filter<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      for (size_t n = 0; n < output_.size(); ++n)
      {
        monitor_->set(static_cast<Variable>(n), [this, n]
                      { return output_[n].read(); });
      }
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Filter<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }
  } // namespace EMT
} // namespace GridKit

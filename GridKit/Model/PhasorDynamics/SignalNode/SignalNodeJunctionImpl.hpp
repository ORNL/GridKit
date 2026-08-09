/**
 * @file SignalNodeJunctionImpl.hpp
 * @brief Implementation of the internal signal-node junction component.
 */
#include <cstddef>
#include <stdexcept>
#include <utility>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeJunction.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    SignalNodeJunction<scalar_type, index_type>::SignalNodeJunction(
        SignalT*                   output,
        RealT                      bias,
        IdxT                       initialization_input_index,
        std::vector<JunctionInput> inputs)
      : output_(output)
    {
      if (output_ == nullptr)
      {
        throw std::invalid_argument("A signal-node junction requires an output signal");
      }

      output_->configureJunction(bias, initialization_input_index, std::move(inputs));
      size_ = 1;
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::allocate()
    {
      if (!allocated_ && output_->linked())
      {
        Log::error() << "SignalNodeJunction: output signal already has a producer\n";
        return 1;
      }

      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      tag_.resize(static_cast<std::size_t>(size_));
      variable_indices_.resize(static_cast<std::size_t>(size_));
      residual_indices_.resize(static_cast<std::size_t>(size_));
      this->allocateExternalVectors(static_cast<IdxT>(output_->junctionInputs().size()));

      this->setVariableIndex(0, 0);
      this->setResidualIndex(0, 0);
      output_->set(y_.getData(), &this->getVariableIndex(0));

      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::verify() const
    {
      int status = 0;

      if (!output_->isJunction())
      {
        Log::error() << "SignalNodeJunction: output is not configured as a junction\n";
        ++status;
      }
      if (!output_->linked())
      {
        Log::error() << "SignalNodeJunction: output signal is not linked\n";
        ++status;
      }

      for (const auto& input : output_->junctionInputs())
      {
        if (!input.signal->linked())
        {
          Log::error() << "SignalNodeJunction: input signal " << input.signal->signalId()
                       << " is not linked\n";
          ++status;
        }
      }

      return status;
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::initialize()
    {
      if (verify() != 0)
      {
        Log::error() << "SignalNodeJunction: cannot initialize with invalid configuration\n";
        return 1;
      }

      output_->initializeJunction();
      yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
      y_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = false;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void SignalNodeJunction<scalar_type, index_type>::gatherExternalVariables()
    {
      const auto& inputs = output_->junctionInputs();
      for (std::size_t i = 0; i < inputs.size(); ++i)
      {
        y_ext_[i]                = inputs[i].signal->read();
        variable_indices_ext_[i] = inputs[i].signal->getVariableIndex();
      }
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::evaluateInternalResidual()
    {
      gatherExternalVariables();

      auto*       f      = f_.getData();
      const auto* y      = y_.getData();
      const auto& inputs = output_->junctionInputs();

      f[0] = y[0] - static_cast<ScalarT>(output_->junctionBias());
      for (std::size_t i = 0; i < inputs.size(); ++i)
      {
        f[0] -= inputs[i].gain * y_ext_[i];
      }

      f_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    int SignalNodeJunction<scalar_type, index_type>::evaluateJacobian()
    {
      gatherExternalVariables();

      if (J_rows_buffer_ == nullptr)
      {
        nnz_ = 1;
        for (const auto index : variable_indices_ext_)
        {
          if (index != INVALID_INDEX<IdxT>)
          {
            ++nnz_;
          }
        }

        J_rows_buffer_ = new IdxT[static_cast<std::size_t>(nnz_)];
        J_cols_buffer_ = new IdxT[static_cast<std::size_t>(nnz_)];
        J_vals_buffer_ = new RealT[static_cast<std::size_t>(nnz_)];

        IdxT entry            = 0;
        J_rows_buffer_[entry] = residual_indices_[0];
        J_cols_buffer_[entry] = variable_indices_[0];
        J_vals_buffer_[entry] = ONE<RealT>;
        ++entry;

        const auto& inputs = output_->junctionInputs();
        for (std::size_t i = 0; i < inputs.size(); ++i)
        {
          if (variable_indices_ext_[i] != INVALID_INDEX<IdxT>)
          {
            J_rows_buffer_[entry] = residual_indices_[0];
            J_cols_buffer_[entry] = variable_indices_ext_[i];
            J_vals_buffer_[entry] = -inputs[i].gain;
            ++entry;
          }
        }
      }

      return this->constructCoo();
    }

  } // namespace PhasorDynamics
} // namespace GridKit

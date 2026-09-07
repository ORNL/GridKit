#pragma once

#include <array>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Bus/BusData.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class KCL : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using Base         = Component<ScalarT, IdxT>;
      using RealT        = typename Base::RealT;
      using SignalT      = typename Base::SignalT;
      using Outputs      = BusOutputs;
      using PhaseSignals = std::array<SignalT*, 3>;
      using PhaseOrder   = std::array<size_t, 3>;

      KCL()
      {
        this->equation_size_ = this->size_ = 3;
        for (auto& signal : voltage_)
          signal.claimProducer();
      }

      SignalT& outputSignal(Outputs output)
      {
        return voltage_.at(static_cast<size_t>(output));
      }

      PhaseSignals voltages(PhaseOrder phases = {0, 1, 2})
      {
        return {&voltage_.at(phases[0]), &voltage_.at(phases[1]), &voltage_.at(phases[2])};
      }

      void attachInput(BusInputs input, SignalT* signal)
      {
        if (this->allocated_)
          throw std::logic_error("KCL inputs cannot change after allocation");
        current_.at(static_cast<size_t>(input)) = signal;
      }

      void assignOutput(Outputs output, SignalT* signal)
      {
        if (this->allocated_)
          throw std::logic_error("KCL outputs cannot change after allocation");
        signal->claimProducer();
        aliases_.at(static_cast<size_t>(output)) = signal;
      }

      IdxT voltagePhase(const SignalT* signal) const
      {
        for (size_t p = 0; p < 3; ++p)
          if (signal == &voltage_[p] || (signal && signal == aliases_[p]))
            return static_cast<IdxT>(p);
        return INVALID_INDEX<IdxT>;
      }

      int setGridKitComponentID(IdxT id) override
      {
        this->gridkit_component_id_ = id;
        return 0;
      }

      int allocate() override
      {
        if (!this->allocated_)
          this->allocateVectors(3);
        this->tag_.resize(3);
        this->variable_indices_.resize(3);
        this->residual_indices_.resize(3);
        this->assignGlobalIndices(0);
        for (IdxT p = 0; p < 3; ++p)
        {
          this->bindSignal(voltage_[static_cast<size_t>(p)], p);
          if (aliases_[static_cast<size_t>(p)])
            this->bindSignal(*aliases_[static_cast<size_t>(p)], p);
        }
        this->allocated_ = true;
        return 0;
      }

      int verify() const override
      {
        int errors = 0;
        for (auto* signal : current_)
          errors += signal && !signal->linked();
        return errors;
      }

      int initializationOrder() const noexcept override
      {
        return 0;
      }

      int initialize(const std::map<Outputs, RealT>& outputs = {})
      {
        this->validateOutputValues(outputs);
        for (size_t p = 0; p < 3; ++p)
        {
          voltage_[p].init(static_cast<ScalarT>(this->outputValue(outputs, static_cast<Outputs>(p), RealT{0})));
          voltage_[p].initDerivative(ScalarT{0});
        }
        this->y_.setDataUpdated();
        this->yp_.setDataUpdated();
        return 0;
      }

      int initializeState(const std::map<std::string, RealT>& values) override
      {
        std::map<Outputs, RealT> outputs;
        for (const auto& [name, value] : values)
        {
          if (name.size() != 2 || name[0] != 'v' || name[1] < 'a' || name[1] > 'c')
            throw std::invalid_argument("Unknown initial bus output: " + name);
          outputs[static_cast<Outputs>(name[1] - 'a')] = value;
        }
        return initialize(outputs);
      }

      int tagDifferentiable() override
      {
        for (size_t p = 0; p < 3; ++p)
          this->tag_[p] = voltage_[p].hasDerivativeCoupling()
                          || (aliases_[p] && aliases_[p]->hasDerivativeCoupling());
        return 0;
      }

      int setAbsoluteTolerance(RealT tolerance) override
      {
        this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        return 0;
      }

      int evaluateInternalResidual() override
      {
        for (size_t p = 0; p < 3; ++p)
          this->f_.getData()[p] = current_[p] ? current_[p]->read() : ScalarT{0};
        this->f_.setDataUpdated();
        return 0;
      }

      int evaluateExternalResidual() override
      {
        return 0;
      }

      int evaluateResidual() override
      {
        return evaluateInternalResidual();
      }

      int evaluateJacobian() override
      {
        std::array<typename SignalT::GradientT, 3> gradients;
        size_t                                     entries = 0;
        for (size_t p = 0; p < 3; ++p)
        {
          if (current_[p])
            current_[p]->appendGradient(gradients[p]);
          entries += gradients[p].size();
        }
        if (entries != capacity_)
        {
          this->resetJacobianStructure();
          delete[] this->J_rows_buffer_;
          delete[] this->J_cols_buffer_;
          delete[] this->J_vals_buffer_;
          this->J_rows_buffer_ = new IdxT[entries];
          this->J_cols_buffer_ = new IdxT[entries];
          this->J_vals_buffer_ = new RealT[entries];
          capacity_            = entries;
        }
        this->nnz_ = 0;
        for (size_t p = 0; p < 3; ++p)
          for (const auto& [column, value] : gradients[p])
          {
            const auto j            = this->nnz_++;
            this->J_rows_buffer_[j] = this->getResidualIndex(static_cast<IdxT>(p));
            this->J_cols_buffer_[j] = column;
            this->J_vals_buffer_[j] = value;
          }
        return entries == 0 ? 0 : this->constructCoo();
      }

    private:
      std::array<SignalT, 3> voltage_;
      PhaseSignals           current_{}, aliases_{};
      size_t                 capacity_{0};
    };
  } // namespace EMT
} // namespace GridKit

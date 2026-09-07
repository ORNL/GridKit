#pragma once

#include <GridKit/Model/EMT/Component/Controller/DCLink/DcLink.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      DcLink<scalar_type, index_type>::DcLink(const ModelDataT& data)
        : capacitance_(parameter<RealT>(data, DcLinkParameters::C)), monitor_(std::make_unique<MonitorT>(data))
      {
        if (capacitance_ <= RealT{0})
          throw std::invalid_argument("DCLink: C must be strictly positive");
        this->equation_size_ = this->size_ = 1;
        voltage_.claimProducer();
        using Mon = typename ModelDataT::MonitorableVariables;
        monitor_->set(Mon::vdc, [this]
                      { return voltage_.read(); });
        monitor_->set(Mon::isrc, [this]
                      { return isrc_->read(); });
        monitor_->set(Mon::idc, [this]
                      { return idc_->read(); });
        monitor_->set(Mon::energy, [this]
                      { const auto v = voltage_.read(); return RealT{0.5} * capacitance_ * v * v; });
      }

      template <typename scalar_type, typename index_type>
      typename DcLink<scalar_type, index_type>::SignalT& DcLink<scalar_type, index_type>::outputSignal(Outputs output)
      {
        if (output != Outputs::vdc)
          throw std::invalid_argument("DCLink: invalid output");
        return voltage_;
      }

      template <typename scalar_type, typename index_type>
      void DcLink<scalar_type, index_type>::attachInput(SignalT* isrc, SignalT* idc)
      {
        if (this->allocated_)
          throw std::logic_error("DCLink inputs cannot change after allocation");
        if (!isrc || !idc)
          throw std::invalid_argument("DCLink requires isrc and idc");
        isrc_ = isrc;
        idc_  = idc;
      }

      template <typename scalar_type, typename index_type>
      void DcLink<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        if (this->allocated_)
          throw std::logic_error("DCLink outputs cannot change after allocation");
        if (output != Outputs::vdc || !signal || alias_)
          throw std::invalid_argument("DCLink: invalid or duplicate output assignment");
        signal->claimProducer();
        alias_ = signal;
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        this->gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::allocate()
      {
        if (!isrc_ || !idc_)
          throw std::invalid_argument("DCLink requires isrc and idc");
        if (!this->allocated_)
          this->allocateVectors(1);
        this->tag_.resize(1);
        this->variable_indices_.resize(1);
        this->residual_indices_.resize(1);
        this->assignGlobalIndices(0);
        this->bindSignal(voltage_, 0);
        if (alias_)
          this->bindSignal(*alias_, 0);
        this->allocateExternalVectors(2, 0);
        this->setExternalVariableSignal(0, isrc_);
        this->setExternalVariableSignal(1, idc_);
        this->allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::verify() const
      {
        return !isrc_ || !idc_ || !isrc_->linked() || !idc_->linked();
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        this->validateOutputValues(outputs);
        voltage_.init(static_cast<ScalarT>(this->outputValue(outputs, Outputs::vdc, RealT{0})));
        voltage_.initDerivative(ScalarT{0});
        this->y_.setDataUpdated();
        this->yp_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void DcLink<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
      {
        this->template parseInitialOutputs<DcLink>(values);
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
      {
        return this->initializeOutputs(*this, values);
      }

      template <typename scalar_type, typename index_type>
      typename Component<scalar_type, index_type>::InitializationPortsT DcLink<scalar_type, index_type>::initializationPorts()
      {
        return {{}, {{"vdc", &voltage_}}, {}};
      }

      template <typename scalar_type, typename index_type>
      void DcLink<scalar_type, index_type>::prepareInitialization(typename Base::InitialStateT& initial)
      {
        const auto outputs = this->template parseInitialOutputs<DcLink>(initial.outputs(*this));
        const auto value   = this->outputValue(outputs, Outputs::vdc, RealT{0});
        initial.provide(voltage_, value);
        if (alias_)
          initial.provide(*alias_, value);
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::evaluateInternalResidual()
      {
        this->f_.getData()[0] = isrc_->read() - idc_->read() - capacitance_ * this->yp_.getData()[0];
        this->f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      int DcLink<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
      {
        typename SignalT::GradientT gradient;
        if (y_scale != ZERO<RealT>)
        {
          isrc_->appendGradient(gradient, y_scale);
          idc_->appendGradient(gradient, -y_scale);
        }
        if (yp_scale != ZERO<RealT>)
          gradient.emplace_back(this->getVariableIndex(0), -capacitance_ * yp_scale);
        const auto entries = gradient.size();
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
        for (const auto& [column, value] : gradient)
        {
          const auto j            = this->nnz_++;
          this->J_rows_buffer_[j] = this->getResidualIndex(0);
          this->J_cols_buffer_[j] = column;
          this->J_vals_buffer_[j] = value;
        }
        return entries == 0 ? 0 : this->constructCoo();
      }

      template <typename scalar_type, typename index_type>
      DcLink<scalar_type, index_type>::~DcLink() = default;

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* DcLink<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit

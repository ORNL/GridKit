#pragma once

#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/ContainerRuntime.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::~Bus() = default;

    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus()
      : kcl_(this->template add<KCLT>("KCL"))
    {
      zero_.bindConstant(ZERO<RealT>);
      for (size_t p = 0; p < 3; ++p)
        this->output(std::string("v") + "abc"[p], outputSignal(static_cast<Outputs>(p)));
    }

    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : Bus()
    {
      monitor_ = std::make_unique<MonitorT>(data);
      for (const auto& [name, Y] : data.shunts)
        addShunt(name, Y);
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    typename Bus<scalar_type, index_type>::NortonT& Bus<scalar_type, index_type>::addNorton(
        std::string name, const YDataT& Y, PhaseSignals incident, RealT scale, PhaseOrder phases)
    {
      for (auto*& signal : incident)
        if (!signal)
          signal = &zero_;
      auto& source = this->template add<NortonT>(name, Y, voltages(phases), scale);
      for (size_t p = 0; p < 3; ++p)
      {
        addCurrent(phases[p], *incident[p]);
        addCurrent(phases[p], source.outputSignal(p), -ONE<RealT>);
        this->input(name + "_inc_" + "abc"[p], *incident[p]);
        this->output(name + "_Ish_" + "abc"[p], source.outputSignal(p));
        shunt_monitors_[phases[p]].push_back(&source.outputSignal(p));
      }
      return source;
    }

    template <typename scalar_type, typename index_type>
    typename Bus<scalar_type, index_type>::NortonT& Bus<scalar_type, index_type>::addShunt(std::string name, const YDataT& Y)
    {
      return addNorton(std::move(name), Y);
    }

    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      int status = kcl_.initialize(outputs);
      this->forEachComponent([&](typename Base::ComponentT& component)
                             {
        if (auto* source = dynamic_cast<NortonT*>(&component); source && status == 0)
          status = source->initialize(); });
      return status;
    }

    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initializeSteadyState(RealT omega)
    {
      int status = 0;
      this->forEachComponent([&](typename Base::ComponentT& component)
                             {
        if (auto* source = dynamic_cast<NortonT*>(&component); source && status == 0)
          status = source->initializeSteadyState(omega); });
      return status;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Bus<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Bus<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      for (size_t p = 0; p < 3; ++p)
      {
        monitor_->set(static_cast<Variable>(p), [this, p]
                      { return outputSignal(static_cast<Outputs>(p)).read(); });
        monitor_->set(static_cast<Variable>(3 + p), [this, p]
                      {
                        ScalarT current{};
                        for (const auto* signal : shunt_monitors_[p])
                          current += signal->read();
                        return current; });
      }
    }
  } // namespace EMT
} // namespace GridKit

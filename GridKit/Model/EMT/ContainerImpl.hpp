#include <algorithm>
#include <stdexcept>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/ComponentLibrary.hpp>
#include <GridKit/Model/EMT/Container.hpp>
#include <GridKit/Model/EMT/ContainerData.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Container<scalar_type, index_type>::Container()
    {
    }

    template <typename scalar_type, typename index_type>
    Container<scalar_type, index_type>::Container(const ModelDataT& data)
    {
      if (!data.inputs.empty())
      {
        throw std::invalid_argument(
            "A standalone Container has no parent scope for its inputs");
      }
      declare(data, {});
      assemble(data);
      validateBoundary();
      refreshLayout();
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::declare(const ModelDataT& data,
                                                     std::string       path)
    {
      path_ = std::move(path);

      // Declare the entire hierarchy before wiring it. This makes sibling
      // outputs available when their peers bind inputs, independent of JSON
      // ordering.
      for (const auto& [name, reference] : data.inputs)
      {
        static_cast<void>(reference);
        declareInput(name);
      }

      for (const auto& signal_data : data.signal)
      {
        auto& signal = addSignal(signal_data.id);
        if (signal_data.value.has_value())
        {
          const auto value = signal_data.value.value();
          signal.setComputed([value]
                             { return static_cast<ScalarT>(value); },
                             [](typename SignalT::GradientT&, RealT) {});
        }
      }

      for (const auto& bus_data : data.bus)
      {
        auto qualified_data = bus_data;
        qualified_data.id   = qualify(bus_data.id);
        add<Bus<ScalarT, IdxT>>(bus_data.id, qualified_data);
      }

      for (const auto& child_data : data.container)
      {
        auto child = std::make_unique<Container>();
        child->declare(child_data, qualify(child_data.id));
        add(child_data.id, std::move(child));
      }

      for (const auto& source_data : data.voltage_source)
      {
        auto qualified_data = source_data;
        qualified_data.id   = qualify(source_data.id);
        add<VoltageSource<ScalarT, IdxT>>(source_data.id, qualified_data);
      }

      for (const auto& source_data : data.dependent_voltage_source)
      {
        auto qualified_data = source_data;
        qualified_data.id   = qualify(source_data.id);
        add<DependentVoltageSource<ScalarT, IdxT>>(source_data.id, qualified_data);
      }

      for (const auto& machine_data : data.machine)
      {
        auto qualified_data = machine_data;
        qualified_data.id   = qualify(machine_data.id);
        add<Machine<ScalarT, IdxT>>(machine_data.id, qualified_data);
      }

      for (const auto& line_data : data.line_lumped)
      {
        auto qualified_data = line_data;
        qualified_data.id   = qualify(line_data.id);
        add<LineLumped<ScalarT, IdxT>>(line_data.id, qualified_data);
      }

      for (const auto& load_data : data.loadz)
      {
        auto qualified_data = load_data;
        qualified_data.id   = qualify(load_data.id);
        add<LoadZ<ScalarT, IdxT>>(load_data.id, qualified_data);
      }

      for (const auto& governor_data : data.gov)
      {
        auto qualified_data = governor_data;
        qualified_data.id   = qualify(governor_data.id);
        add<Controller::Tgov1<ScalarT, IdxT>>(governor_data.id, qualified_data);
      }

      for (const auto& exciter_data : data.exciter)
      {
        auto qualified_data = exciter_data;
        qualified_data.id   = qualify(exciter_data.id);
        add<Controller::Ieeet1<ScalarT, IdxT>>(exciter_data.id, qualified_data);
      }

      for (const auto& switch_data : data.sw)
      {
        auto qualified_data = switch_data;
        qualified_data.id   = qualify(switch_data.id);
        add<Switch<ScalarT, IdxT>>(switch_data.id, qualified_data);
      }

      for (const auto& model_data : data.pwm)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Controller::Pwm<ScalarT, IdxT>>(model_data.id, qualified_data);
      }

      for (const auto& model_data : data.converter)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Converter<ScalarT, IdxT>>(model_data.id, qualified_data);
      }

      // Outputs are aliases of already-declared internal endpoints. Publishing
      // them bottom-up lets a parent bind one child's input to a sibling's
      // output before any leaf component is wired.
      for (const auto& [name, reference] : data.outputs)
      {
        const auto value = resolveOutput(reference);
        if (std::holds_alternative<SignalT*>(value))
        {
          output(name, *std::get<SignalT*>(value));
        }
        else
        {
          output(name, *std::get<Port3T*>(value));
        }
      }

      refreshLayout();
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::assemble(const ModelDataT& data)
    {
      // A child input is a local name inside the child and a reference in this
      // scope. Bind every child boundary before wiring anything below it.
      for (const auto& child_data : data.container)
      {
        auto& child = component<Container>(child_data.id);
        for (const auto& [name, reference] : child_data.inputs)
        {
          child.bindInput(name, endpoint(reference));
        }
      }

      for (const auto& child_data : data.container)
      {
        component<Container>(child_data.id).assemble(child_data);
      }

      wire(data);

      for (const auto& child_data : data.container)
      {
        component<Container>(child_data.id).validateBoundary();
      }
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::wire(const ModelDataT& data)
    {
      for (const auto& model_data : data.pwm)
      {
        auto& model = component<Controller::Pwm<ScalarT, IdxT>>(model_data.id);
        if (model_data.outputs.contains(Controller::PwmOutputs::sa))
        {
          model.assignOutput(0, &signal(model_data.outputs.at(Controller::PwmOutputs::sa)));
        }
        if (model_data.outputs.contains(Controller::PwmOutputs::sb))
        {
          model.assignOutput(1, &signal(model_data.outputs.at(Controller::PwmOutputs::sb)));
        }
        if (model_data.outputs.contains(Controller::PwmOutputs::sc))
        {
          model.assignOutput(2, &signal(model_data.outputs.at(Controller::PwmOutputs::sc)));
        }
      }

      for (const auto& model_data : data.converter)
      {
        auto& model = component<Converter<ScalarT, IdxT>>(model_data.id);
        model.attachInput(&source(model_data.inputs.at(ConverterInputs::sa)),
                          &source(model_data.inputs.at(ConverterInputs::sb)),
                          &source(model_data.inputs.at(ConverterInputs::sc)),
                          &source(model_data.inputs.at(ConverterInputs::vdc)));
        if (model_data.outputs.contains(ConverterOutputs::voa))
        {
          model.assignOutput(0, &signal(model_data.outputs.at(ConverterOutputs::voa)));
        }
        if (model_data.outputs.contains(ConverterOutputs::vob))
        {
          model.assignOutput(1, &signal(model_data.outputs.at(ConverterOutputs::vob)));
        }
        if (model_data.outputs.contains(ConverterOutputs::voc))
        {
          model.assignOutput(2, &signal(model_data.outputs.at(ConverterOutputs::voc)));
        }
      }

      for (const auto& source_data : data.voltage_source)
      {
        auto&       source_model = component<VoltageSource<ScalarT, IdxT>>(source_data.id);
        const auto& bus_ref      = source_data.inputs.at(VoltageSourceInputs::bus);
        source_model.getSignals().template attachPort<VoltageSourceExternalVariables::VA>(
            &port(bus_ref));
      }

      for (const auto& source_data : data.dependent_voltage_source)
      {
        auto&       source_model = component<DependentVoltageSource<ScalarT, IdxT>>(source_data.id);
        const auto& bus_ref      = source_data.inputs.at(DependentVoltageSourceInputs::bus);
        source_model.getSignals().template attachPort<DependentVoltageSourceExternalVariables::VA>(
            &port(bus_ref));

        if (source_data.inputs.contains(DependentVoltageSourceInputs::ea))
        {
          source_model.getSignals().template attachSignal<DependentVoltageSourceExternalVariables::EA>(
              &source(source_data.inputs.at(DependentVoltageSourceInputs::ea)));
        }
        if (source_data.inputs.contains(DependentVoltageSourceInputs::eb))
        {
          source_model.getSignals().template attachSignal<DependentVoltageSourceExternalVariables::EB>(
              &source(source_data.inputs.at(DependentVoltageSourceInputs::eb)));
        }
        if (source_data.inputs.contains(DependentVoltageSourceInputs::ec))
        {
          source_model.getSignals().template attachSignal<DependentVoltageSourceExternalVariables::EC>(
              &source(source_data.inputs.at(DependentVoltageSourceInputs::ec)));
        }
      }

      for (const auto& machine_data : data.machine)
      {
        auto& machine_model = component<Machine<ScalarT, IdxT>>(machine_data.id);
        machine_model.getSignals().template attachPort<MachineExternalVariables::VA>(
            &port(machine_data.inputs.at(MachineInputs::bus)));

        if (machine_data.outputs.contains(MachineOutputs::speed))
        {
          machine_model.getSignals().template assignSignal<MachineInternalVariables::OMEGA>(
              &signal(machine_data.outputs.at(MachineOutputs::speed)));
        }
        if (machine_data.inputs.contains(MachineInputs::pm))
        {
          machine_model.getSignals().template attachSignal<MachineExternalVariables::PM>(
              &source(machine_data.inputs.at(MachineInputs::pm)));
        }
        if (machine_data.inputs.contains(MachineInputs::efd))
        {
          machine_model.getSignals().template attachSignal<MachineExternalVariables::EFD>(
              &source(machine_data.inputs.at(MachineInputs::efd)));
        }
      }

      for (const auto& line_data : data.line_lumped)
      {
        auto& line_model = component<LineLumped<ScalarT, IdxT>>(line_data.id);
        line_model.getSignals().template attachPort<LineLumpedExternalVariables::V1A>(
            &port(line_data.inputs.at(LineLumpedInputs::bus1)));
        line_model.getSignals().template attachPort<LineLumpedExternalVariables::V2A>(
            &port(line_data.inputs.at(LineLumpedInputs::bus2)));
      }

      for (const auto& load_data : data.loadz)
      {
        auto& load_model = component<LoadZ<ScalarT, IdxT>>(load_data.id);
        load_model.getSignals().template attachPort<LoadZExternalVariables::VA>(
            &port(load_data.inputs.at(LoadZInputs::bus)));
      }

      for (const auto& governor_data : data.gov)
      {
        auto& governor_model = component<Controller::Tgov1<ScalarT, IdxT>>(governor_data.id);
        if (governor_data.inputs.contains(Controller::Tgov1Inputs::speed))
        {
          governor_model.getSignals().template attachSignal<Controller::Tgov1ExternalVariables::OMEGA>(
              &source(governor_data.inputs.at(Controller::Tgov1Inputs::speed)));
        }
        if (governor_data.inputs.contains(Controller::Tgov1Inputs::pref))
        {
          governor_model.getSignals().template attachSignal<Controller::Tgov1ExternalVariables::PREF>(
              &source(governor_data.inputs.at(Controller::Tgov1Inputs::pref)));
        }
        if (governor_data.outputs.contains(Controller::Tgov1Outputs::pmech))
        {
          governor_model.getSignals().template assignSignal<Controller::Tgov1InternalVariables::PM>(
              &signal(governor_data.outputs.at(Controller::Tgov1Outputs::pmech)));
        }
      }

      for (const auto& exciter_data : data.exciter)
      {
        using Inputs        = Controller::Ieeet1Inputs;
        using External      = Controller::Ieeet1ExternalVariables;
        using Internal      = Controller::Ieeet1InternalVariables;
        auto& exciter_model = component<Controller::Ieeet1<ScalarT, IdxT>>(exciter_data.id);
        auto& signals       = exciter_model.getSignals();
        signals.template attachPort<External::VA>(&port(exciter_data.inputs.at(Inputs::bus)));
        if (exciter_data.inputs.contains(Inputs::speed))
        {
          signals.template attachSignal<External::OMEGA>(&source(exciter_data.inputs.at(Inputs::speed)));
        }
        if (exciter_data.inputs.contains(Inputs::vref))
        {
          signals.template attachSignal<External::VREF>(&source(exciter_data.inputs.at(Inputs::vref)));
        }
        if (exciter_data.inputs.contains(Inputs::vs))
        {
          signals.template attachSignal<External::VS>(&source(exciter_data.inputs.at(Inputs::vs)));
        }
        if (exciter_data.inputs.contains(Inputs::vuel))
        {
          signals.template attachSignal<External::VUEL>(&source(exciter_data.inputs.at(Inputs::vuel)));
        }
        if (exciter_data.inputs.contains(Inputs::voel))
        {
          signals.template attachSignal<External::VOEL>(&source(exciter_data.inputs.at(Inputs::voel)));
        }
        if (exciter_data.outputs.contains(Controller::Ieeet1Outputs::efd))
        {
          signals.template assignSignal<Internal::EFD>(
              &signal(exciter_data.outputs.at(Controller::Ieeet1Outputs::efd)));
        }
      }

      for (const auto& switch_data : data.sw)
      {
        auto& switch_model = component<Switch<ScalarT, IdxT>>(switch_data.id);
        switch_model.getSignals().template attachPort<SwitchExternalVariables::V1A>(
            &port(switch_data.inputs.at(SwitchInputs::bus1)));
        switch_model.getSignals().template attachPort<SwitchExternalVariables::V2A>(
            &port(switch_data.inputs.at(SwitchInputs::bus2)));
      }
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::validateName(std::string_view name,
                                                          std::string_view kind)
    {
      if (name.empty())
      {
        throw std::invalid_argument(std::string(kind) + " name must not be empty");
      }
      if (name.find('.') != std::string_view::npos)
      {
        throw std::invalid_argument(std::string(kind) + " name \"" + std::string(name)
                                    + "\" must not contain '.'");
      }
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::ComponentT&
    Container<scalar_type, index_type>::add(std::string                 id,
                                            std::unique_ptr<ComponentT> child)
    {
      if (allocated_ || bound_)
      {
        throw std::logic_error("A Container cannot change after allocation");
      }
      validateName(id, "Component");
      if (child == nullptr)
      {
        throw std::invalid_argument("A Container cannot own a null child");
      }
      if (children_by_id_.contains(id) || signals_by_id_.contains(id)
          || input_names_.contains(id))
      {
        throw std::invalid_argument("Duplicate local name \"" + id + "\"");
      }

      auto* ptr = child.get();
      ptr->setGridKitComponentID(static_cast<IdxT>(children_.size()));
      children_by_id_.emplace(id, ptr);
      children_.push_back(std::move(child));
      refreshLayout();
      return *ptr;
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::addSignal(std::string id)
    {
      if (allocated_ || bound_)
      {
        throw std::logic_error("A Container cannot change after allocation");
      }
      validateName(id, "Signal");
      if (signals_by_id_.contains(id) || children_by_id_.contains(id)
          || input_names_.contains(id))
      {
        throw std::invalid_argument("Duplicate local name \"" + id + "\"");
      }

      auto  value = std::make_unique<SignalT>(id);
      auto* ptr   = value.get();
      signals_by_id_.emplace(id, ptr);
      signals_.push_back(std::move(value));
      return *ptr;
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::ComponentT&
    Container<scalar_type, index_type>::component(std::string_view path)
    {
      const auto dot = path.find('.');
      if (dot == std::string_view::npos)
      {
        const auto found = children_by_id_.find(path);
        if (found == children_by_id_.end())
        {
          throw std::invalid_argument("Unknown component \"" + std::string(path) + "\"");
        }
        return *found->second;
      }
      return childContainer(path.substr(0, dot)).component(path.substr(dot + 1));
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::ComponentT&
    Container<scalar_type, index_type>::component(std::string_view path) const
    {
      const auto dot = path.find('.');
      if (dot == std::string_view::npos)
      {
        const auto found = children_by_id_.find(path);
        if (found == children_by_id_.end())
        {
          throw std::invalid_argument("Unknown component \"" + std::string(path) + "\"");
        }
        return *found->second;
      }
      return childContainer(path.substr(0, dot)).component(path.substr(dot + 1));
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::ComponentT&
    Container<scalar_type, index_type>::component(IdxT local_index)
    {
      return *children_.at(static_cast<size_t>(local_index));
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::ComponentT&
    Container<scalar_type, index_type>::component(IdxT local_index) const
    {
      return *children_.at(static_cast<size_t>(local_index));
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::signal(std::string_view path)
    {
      const auto dot = path.find('.');
      if (dot != std::string_view::npos)
      {
        return childContainer(path.substr(0, dot)).signal(path.substr(dot + 1));
      }
      const auto found = signals_by_id_.find(path);
      if (found == signals_by_id_.end())
      {
        throw std::invalid_argument("Unknown signal \"" + std::string(path) + "\"");
      }
      return *found->second;
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::signal(std::string_view path) const
    {
      const auto dot = path.find('.');
      if (dot != std::string_view::npos)
      {
        return childContainer(path.substr(0, dot)).signal(path.substr(dot + 1));
      }
      const auto found = signals_by_id_.find(path);
      if (found == signals_by_id_.end())
      {
        throw std::invalid_argument("Unknown signal \"" + std::string(path) + "\"");
      }
      return *found->second;
    }

    template <typename scalar_type, typename index_type>
    Container<scalar_type, index_type>&
    Container<scalar_type, index_type>::childContainer(std::string_view id)
    {
      auto* child = dynamic_cast<Container*>(&component(id));
      if (child == nullptr)
      {
        throw std::invalid_argument("Component \"" + std::string(id) + "\" is not a Container");
      }
      return *child;
    }

    template <typename scalar_type, typename index_type>
    const Container<scalar_type, index_type>&
    Container<scalar_type, index_type>::childContainer(std::string_view id) const
    {
      auto* child = dynamic_cast<const Container*>(&component(id));
      if (child == nullptr)
      {
        throw std::invalid_argument("Component \"" + std::string(id) + "\" is not a Container");
      }
      return *child;
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::input(std::string name, SignalT& signal_value)
    {
      if (allocated_ || bound_)
      {
        throw std::logic_error("A Container cannot change after allocation");
      }
      validateName(name, "Input");
      if (!input_names_.contains(name))
      {
        declareInput(name);
      }
      bindInput(name, &signal_value);
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::input(std::string name, Port3T& port_value)
    {
      if (allocated_ || bound_)
      {
        throw std::logic_error("A Container cannot change after allocation");
      }
      validateName(name, "Input");
      if (!input_names_.contains(name))
      {
        declareInput(name);
      }
      bindInput(name, &port_value);
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::declareInput(std::string name)
    {
      validateName(name, "Input");
      if (input_names_.contains(name) || outputs_.contains(name)
          || signals_by_id_.contains(name) || children_by_id_.contains(name))
      {
        throw std::invalid_argument("Duplicate local or boundary name \"" + name + "\"");
      }
      input_names_.insert(std::move(name));
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::bindInput(std::string_view name,
                                                       Endpoint         value)
    {
      if (!input_names_.contains(name))
      {
        throw std::invalid_argument("Unknown Container input \"" + std::string(name) + "\"");
      }
      if (!inputs_.emplace(std::string(name), value).second)
      {
        throw std::invalid_argument("Container input \"" + std::string(name)
                                    + "\" is already bound");
      }
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::Endpoint
    Container<scalar_type, index_type>::inputEndpoint(std::string_view name) const
    {
      const auto found = inputs_.find(name);
      if (found == inputs_.end())
      {
        throw std::invalid_argument("Unbound Container input \"" + std::string(name) + "\"");
      }
      return found->second;
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::inputSignal(std::string_view name)
    {
      const auto value = inputEndpoint(name);
      if (!std::holds_alternative<SignalT*>(value))
      {
        throw std::invalid_argument("Container input \"" + std::string(name)
                                    + "\" is not a scalar signal");
      }
      return *std::get<SignalT*>(value);
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::inputSignal(std::string_view name) const
    {
      const auto value = inputEndpoint(name);
      if (!std::holds_alternative<SignalT*>(value))
      {
        throw std::invalid_argument("Container input \"" + std::string(name)
                                    + "\" is not a scalar signal");
      }
      return *std::get<SignalT*>(value);
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::Port3T&
    Container<scalar_type, index_type>::inputPort(std::string_view name)
    {
      const auto value = inputEndpoint(name);
      if (!std::holds_alternative<Port3T*>(value))
      {
        throw std::invalid_argument("Container input \"" + std::string(name)
                                    + "\" is not a three-phase port");
      }
      return *std::get<Port3T*>(value);
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::Port3T&
    Container<scalar_type, index_type>::inputPort(std::string_view name) const
    {
      const auto value = inputEndpoint(name);
      if (!std::holds_alternative<Port3T*>(value))
      {
        throw std::invalid_argument("Container input \"" + std::string(name)
                                    + "\" is not a three-phase port");
      }
      return *std::get<Port3T*>(value);
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::output(std::string name, SignalT& signal_value)
    {
      if (allocated_ || bound_)
      {
        throw std::logic_error("A Container cannot change after allocation");
      }
      validateName(name, "Output");
      if (input_names_.contains(name) || outputs_.contains(name))
      {
        throw std::invalid_argument("Duplicate boundary name \"" + name + "\"");
      }
      outputs_.emplace(std::move(name), &signal_value);
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::output(std::string name, Port3T& port_value)
    {
      if (allocated_ || bound_)
      {
        throw std::logic_error("A Container cannot change after allocation");
      }
      validateName(name, "Output");
      if (input_names_.contains(name) || outputs_.contains(name))
      {
        throw std::invalid_argument("Duplicate boundary name \"" + name + "\"");
      }
      outputs_.emplace(std::move(name), &port_value);
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::outputSignal(std::string_view name)
    {
      const auto found = outputs_.find(name);
      if (found == outputs_.end() || !std::holds_alternative<SignalT*>(found->second))
      {
        throw std::invalid_argument("Container output \"" + std::string(name)
                                    + "\" is not a scalar signal");
      }
      return *std::get<SignalT*>(found->second);
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::outputSignal(std::string_view name) const
    {
      const auto found = outputs_.find(name);
      if (found == outputs_.end() || !std::holds_alternative<SignalT*>(found->second))
      {
        throw std::invalid_argument("Container output \"" + std::string(name)
                                    + "\" is not a scalar signal");
      }
      return *std::get<SignalT*>(found->second);
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::Port3T&
    Container<scalar_type, index_type>::outputPort(std::string_view name)
    {
      const auto found = outputs_.find(name);
      if (found == outputs_.end() || !std::holds_alternative<Port3T*>(found->second))
      {
        throw std::invalid_argument("Container output \"" + std::string(name)
                                    + "\" is not a three-phase port");
      }
      return *std::get<Port3T*>(found->second);
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::Port3T&
    Container<scalar_type, index_type>::outputPort(std::string_view name) const
    {
      const auto found = outputs_.find(name);
      if (found == outputs_.end() || !std::holds_alternative<Port3T*>(found->second))
      {
        throw std::invalid_argument("Container output \"" + std::string(name)
                                    + "\" is not a three-phase port");
      }
      return *std::get<Port3T*>(found->second);
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::source(std::string_view reference)
    {
      const auto value = endpoint(reference);
      if (!std::holds_alternative<SignalT*>(value))
      {
        throw std::invalid_argument("Endpoint \"" + std::string(reference)
                                    + "\" is not a scalar signal");
      }
      return *std::get<SignalT*>(value);
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::Port3T&
    Container<scalar_type, index_type>::port(std::string_view reference)
    {
      const auto value = endpoint(reference);
      if (!std::holds_alternative<Port3T*>(value))
      {
        throw std::invalid_argument("Endpoint \"" + std::string(reference)
                                    + "\" is not a three-phase port");
      }
      return *std::get<Port3T*>(value);
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::Endpoint
    Container<scalar_type, index_type>::endpoint(std::string_view reference)
    {
      const auto dot = reference.find('.');
      if (dot != std::string_view::npos)
      {
        return childContainer(reference.substr(0, dot)).outputEndpoint(reference.substr(dot + 1));
      }
      if (const auto found = inputs_.find(reference); found != inputs_.end())
      {
        return found->second;
      }
      if (const auto found = signals_by_id_.find(reference); found != signals_by_id_.end())
      {
        return found->second;
      }
      if (const auto found = children_by_id_.find(reference); found != children_by_id_.end())
      {
        if (auto* bus = dynamic_cast<Bus<ScalarT, IdxT>*>(found->second); bus != nullptr)
        {
          return &bus->voltagePort();
        }
      }
      throw std::invalid_argument("Unknown endpoint \"" + std::string(reference) + "\"");
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::Endpoint
    Container<scalar_type, index_type>::resolveOutput(std::string_view reference)
    {
      const auto dot = reference.find('.');
      if (dot != std::string_view::npos)
      {
        return childContainer(reference.substr(0, dot)).outputEndpoint(reference.substr(dot + 1));
      }

      if (const auto found = signals_by_id_.find(reference); found != signals_by_id_.end())
      {
        return found->second;
      }
      if (const auto found = children_by_id_.find(reference); found != children_by_id_.end())
      {
        auto* bus = dynamic_cast<Bus<ScalarT, IdxT>*>(found->second);
        if (bus != nullptr)
        {
          return &bus->voltagePort();
        }
      }
      throw std::invalid_argument("Output reference \"" + std::string(reference)
                                  + "\" is not a local signal, Bus, or child output");
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::Endpoint
    Container<scalar_type, index_type>::outputEndpoint(std::string_view name)
    {
      const auto found = outputs_.find(name);
      if (found == outputs_.end())
      {
        throw std::invalid_argument("Unknown Container output \"" + std::string(name) + "\"");
      }
      return found->second;
    }

    template <typename scalar_type, typename index_type>
    std::string Container<scalar_type, index_type>::qualify(std::string_view local_name) const
    {
      if (path_.empty())
      {
        return std::string(local_name);
      }
      return path_ + "." + std::string(local_name);
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::validateBoundary() const
    {
      for (const auto& name : input_names_)
      {
        if (!inputs_.contains(name))
        {
          throw std::invalid_argument("Container input \"" + name + "\" is not bound");
        }
      }

      for (const auto& [name, value] : outputs_)
      {
        if (std::holds_alternative<SignalT*>(value)
            && !std::get<SignalT*>(value)->hasProducer())
        {
          throw std::invalid_argument("Container output \"" + name
                                      + "\" has no internal producer");
        }
      }
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::refreshLayout()
    {
      offsets_.clear();
      size_ = 0;
      for (auto& child : children_)
      {
        offsets_.push_back(size_);
        size_ += child->size();
      }
    }

    template <typename scalar_type, typename index_type>
    index_type Container<scalar_type, index_type>::size()
    {
      refreshLayout();
      return size_;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::bind(VectorT& y,
                                                 VectorT& yp,
                                                 VectorT& f,
                                                 VectorT& abs_tol,
                                                 IdxT     offset)
    {
      refreshLayout();
      if (y.getSize() < offset + size_ || yp.getSize() < offset + size_
          || f.getSize() < offset + size_ || abs_tol.getSize() < offset + size_)
      {
        Log::error() << "Container::bind - system vectors are smaller than offset + size = "
                     << offset + size_ << '\n';
        return 1;
      }

      auto* y_data       = y.getData(memory::HOST);
      auto* yp_data      = yp.getData(memory::HOST);
      auto* f_data       = f.getData(memory::HOST);
      auto* abs_tol_data = abs_tol.getData(memory::HOST);
      if (size_ != 0 && (y_data == nullptr || yp_data == nullptr || f_data == nullptr || abs_tol_data == nullptr))
      {
        Log::error() << "Container::bind - system vector data is null or stale\n";
        return 1;
      }

      const int y_status = y_.setData(y_data == nullptr ? nullptr : y_data + offset,
                                      size_,
                                      memory::HOST);
      const int yp_status = yp_.setData(yp_data == nullptr ? nullptr : yp_data + offset,
                                        size_,
                                        memory::HOST);
      const int f_status = f_.setData(f_data == nullptr ? nullptr : f_data + offset,
                                      size_,
                                      memory::HOST);
      const int abs_tol_status = abs_tol_.setData(
          abs_tol_data == nullptr ? nullptr : abs_tol_data + offset, size_, memory::HOST);
      if (y_status != 0 || yp_status != 0 || f_status != 0 || abs_tol_status != 0)
      {
        Log::error() << "Container::bind - failed to bind vectors to system storage\n";
        return 1;
      }

      bound_     = true;
      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::allocate()
    {
      validateBoundary();
      refreshLayout();
      if (!bound_ && !allocated_)
      {
        this->allocateVectors(size_);
      }

      if (y_.getSize() != size_ || yp_.getSize() != size_
          || f_.getSize() != size_ || abs_tol_.getSize() != size_)
      {
        throw std::runtime_error("Container vector sizes do not match its child layout");
      }

      tag_.resize(static_cast<size_t>(size_));
      variable_indices_.resize(static_cast<size_t>(size_));
      residual_indices_.resize(static_cast<size_t>(size_));

      for (size_t i = 0; i < children_.size(); ++i)
      {
        auto& child = children_[i];
        if (child->bind(y_, yp_, f_, abs_tol_, offsets_[i]) != 0
            || child->allocate() != 0)
        {
          throw std::runtime_error("Failed to allocate a Container child");
        }
      }

      assignGlobalIndices(0);
      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::assignGlobalIndices(IdxT first)
    {
      for (IdxT j = 0; j < size_; ++j)
      {
        variable_indices_[static_cast<size_t>(j)] = first + j;
        residual_indices_[static_cast<size_t>(j)] = first + j;
      }
      for (size_t i = 0; i < children_.size(); ++i)
      {
        children_[i]->assignGlobalIndices(first + offsets_[i]);
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::verify() const
    {
      int errors = 0;
      for (const auto& child : children_)
      {
        errors += child->verify();
      }
      return errors;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::initialize()
    {
      std::vector<ComponentT*> leaves;
      forEachComponent([&leaves](ComponentT& component)
                       {
        if (dynamic_cast<Container*>(&component) == nullptr)
        {
          leaves.push_back(&component);
        } });
      std::stable_sort(leaves.begin(), leaves.end(), [](const auto* lhs, const auto* rhs)
                       { return lhs->initializationOrder() < rhs->initializationOrder(); });

      int status = 0;
      for (auto* leaf : leaves)
      {
        status += leaf->initialize();
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return status;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::tagDifferentiable()
    {
      std::fill(tag_.begin(), tag_.end(), false);
      for (size_t i = 0; i < children_.size(); ++i)
      {
        auto& child = children_[i];
        child->tagDifferentiable();
        const auto& child_tag = child->tag();
        for (size_t j = 0; j < child_tag.size(); ++j)
        {
          tag_[static_cast<size_t>(offsets_[i]) + j] = child_tag[j];
        }
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      int status = 0;
      for (auto& child : children_)
      {
        status += child->setAbsoluteTolerance(rel_tol);
      }
      abs_tol_.setDataUpdated();
      return status;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::evaluateInternalResidual()
    {
      for (auto& child : children_)
      {
        const int status = child->evaluateInternalResidual();
        if (status != 0)
        {
          return status;
        }
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::evaluateExternalResidual()
    {
      for (auto& child : children_)
      {
        const int status = child->evaluateExternalResidual();
        if (status != 0)
        {
          return status;
        }
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::evaluateResidual()
    {
      const int internal_status = evaluateInternalResidual();
      if (internal_status != 0)
      {
        return internal_status;
      }
      const int external_status = evaluateExternalResidual();
      f_.setDataUpdated();
      return external_status;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::evaluateJacobian()
    {
      size_t required = 0;
      for (auto& child : children_)
      {
        const int status = child->evaluateJacobian();
        if (status != 0)
        {
          return status;
        }
        if (auto* jacobian = child->getCooJacobian(); jacobian != nullptr)
        {
          required += static_cast<size_t>(jacobian->getNnz());
        }
      }

      if (required == 0)
      {
        delete coo_jac_;
        coo_jac_ = nullptr;
        nnz_     = 0;
        return 0;
      }

      if (required > jacobian_capacity_)
      {
        delete coo_jac_;
        coo_jac_ = nullptr;
        delete[] J_rows_buffer_;
        delete[] J_cols_buffer_;
        delete[] J_vals_buffer_;
        J_rows_buffer_     = new IdxT[required];
        J_cols_buffer_     = new IdxT[required];
        J_vals_buffer_     = new RealT[required];
        jacobian_capacity_ = required;
      }

      nnz_ = 0;
      for (const auto& child : children_)
      {
        auto* jacobian = child->getCooJacobian();
        if (jacobian == nullptr)
        {
          continue;
        }
        for (IdxT j = 0; j < jacobian->getNnz(); ++j)
        {
          J_rows_buffer_[nnz_] = jacobian->getRowData()[j];
          J_cols_buffer_[nnz_] = jacobian->getColData()[j];
          J_vals_buffer_[nnz_] = jacobian->getValues()[j];
          ++nnz_;
        }
      }

      this->constructCoo();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    bool Container<scalar_type, index_type>::hasJacobian()
    {
      for (auto& child : children_)
      {
        if (!child->hasJacobian())
        {
          return false;
        }
      }
      return true;
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::updateTime(RealT t, RealT a)
    {
      time_  = t;
      alpha_ = a;
      for (auto& child : children_)
      {
        child->updateTime(t, a);
      }
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::resetJacobianStructure()
    {
      ComponentT::resetJacobianStructure();
      for (auto& child : children_)
      {
        child->resetJacobianStructure();
      }
    }
  } // namespace EMT
} // namespace GridKit

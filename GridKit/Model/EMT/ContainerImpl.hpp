#include <algorithm>
#include <stdexcept>
#include <tuple>

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/ComponentLibrary.hpp>
#include <GridKit/Model/EMT/ContainerData.hpp>
#include <GridKit/Model/EMT/ContainerRuntime.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Container<scalar_type, index_type>::Container(const ModelDataT& data)
    {
      if (!data.inputs.empty())
      {
        throw std::invalid_argument(
            "A standalone Container has no parent scope for its inputs");
      }
      declare(data, {});
      std::vector<BusT*> buses;
      forEachComponent([&](ComponentT& component)
                       {
                         if (auto* bus = dynamic_cast<BusT*>(&component))
                           buses.push_back(bus); });
      assemble(data, buses);
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
          signal.bindConstant(*signal_data.value);
        }
      }

      for (const auto& bus_data : data.bus)
      {
        auto qualified_data = bus_data;
        qualified_data.id   = qualify(bus_data.id);
        auto& bus           = add<Bus<ScalarT, IdxT>>(bus_data.id, qualified_data);
        for (const auto& [output, reference] : bus_data.outputs)
          bus.assignOutput(output, &signal(reference));
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

      for (const auto& filter_data : data.filter)
      {
        auto qualified_data = filter_data;
        qualified_data.id   = qualify(filter_data.id);
        add<Filter<ScalarT, IdxT>>(filter_data.id, qualified_data);
      }

      for (const auto& machine_data : data.machine)
      {
        auto qualified_data = machine_data;
        qualified_data.id   = qualify(machine_data.id);
        add<Machine<ScalarT, IdxT>>(machine_data.id, qualified_data);
      }

      for (const auto& source_data : data.regfma)
      {
        auto qualified_data = source_data;
        qualified_data.id   = qualify(source_data.id);
        add<Regfma<ScalarT, IdxT>>(source_data.id, qualified_data);
      }

      for (const auto& line_data : data.line_lumped)
      {
        auto qualified_data = line_data;
        qualified_data.id   = qualify(line_data.id);
        add<LineLumped<ScalarT, IdxT>>(line_data.id, qualified_data);
      }

      for (const auto& line_data : data.line_distributed)
      {
        auto qualified_data = line_data;
        qualified_data.id   = qualify(line_data.id);
        add<LineDistributed<ScalarT, IdxT>>(line_data.id, qualified_data);
      }

      for (const auto& load_data : data.loadz)
      {
        auto qualified_data = load_data;
        qualified_data.id   = qualify(load_data.id);
        add<LoadZ<ScalarT, IdxT>>(load_data.id, qualified_data);
      }

      for (const auto& stabilizer_data : data.ieeest)
      {
        auto qualified_data = stabilizer_data;
        qualified_data.id   = qualify(stabilizer_data.id);
        add<Controller::Ieeest<ScalarT, IdxT>>(stabilizer_data.id, qualified_data);
      }

      for (const auto& governor_data : data.gastpti)
      {
        auto qualified_data = governor_data;
        qualified_data.id   = qualify(governor_data.id);
        add<Controller::GastPti<ScalarT, IdxT>>(governor_data.id, qualified_data);
      }

      for (const auto& governor_data : data.gov)
      {
        auto qualified_data = governor_data;
        qualified_data.id   = qualify(governor_data.id);
        add<Controller::Tgov1<ScalarT, IdxT>>(governor_data.id, qualified_data);
      }

      for (const auto& exciter_data : data.sexs_pti)
      {
        auto qualified_data = exciter_data;
        qualified_data.id   = qualify(exciter_data.id);
        add<Controller::SexsPti<ScalarT, IdxT>>(exciter_data.id, qualified_data);
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

      for (const auto& transformer_data : data.transformer)
      {
        auto qualified_data = transformer_data;
        qualified_data.id   = qualify(transformer_data.id);
        add<Transformer<ScalarT, IdxT>>(transformer_data.id, qualified_data);
      }

      for (const auto& model_data : data.inner_current_control)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Controller::InnerCurrentControl<ScalarT, IdxT>>(model_data.id, qualified_data);
      }

      for (const auto& model_data : data.outer_voltage_control)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Controller::OuterVoltageControl<ScalarT, IdxT>>(model_data.id, qualified_data);
      }
      for (const auto& model_data : data.outer_power_control)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Controller::OuterPowerControl<ScalarT, IdxT>>(model_data.id, qualified_data);
      }

      for (const auto& model_data : data.park)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Park<ScalarT, IdxT>>(model_data.id, qualified_data);
      }

      for (const auto& model_data : data.pll)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Pll<ScalarT, IdxT>>(model_data.id, qualified_data);
      }

      for (const auto& model_data : data.pwm)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Controller::Pwm<ScalarT, IdxT>>(model_data.id, qualified_data);
      }

      for (const auto& model_data : data.dc_link)
      {
        auto qualified_data = model_data;
        qualified_data.id   = qualify(model_data.id);
        add<Controller::DcLink<ScalarT, IdxT>>(model_data.id, qualified_data);
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
        output(name, *value);
      }

      refreshLayout();
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::assemble(const ModelDataT& data, const std::vector<BusT*>& buses)
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
        component<Container>(child_data.id).assemble(child_data, buses);
      }

      wire(data, buses);

      for (const auto& child_data : data.container)
      {
        component<Container>(child_data.id).validateBoundary();
      }
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::wire(const ModelDataT& data, const std::vector<BusT*>& buses)
    {
      auto terminal = [&](const auto& inputs, auto first)
      {
        std::array<SignalT*, 3> voltage;
        for (size_t p = 0; p < 3; ++p)
          voltage[p] = &source(inputs.at(static_cast<decltype(first)>(static_cast<size_t>(first) + p)));
        for (auto* bus : buses)
        {
          typename BusT::PhaseOrder phases;
          bool                      found = true;
          for (size_t p = 0; p < 3; ++p)
          {
            const IdxT phase = bus->voltagePhase(voltage[p]);
            found            = found && phase != INVALID_INDEX<IdxT>;
            phases[p]        = static_cast<size_t>(phase);
          }
          if (found)
            return std::make_tuple(bus, phases, voltage);
        }
        throw std::invalid_argument("Terminal phases must belong to one Bus");
      };

      for (const auto& bus_data : data.bus)
      {
        auto& bus = component<Bus<ScalarT, IdxT>>(bus_data.id);
        for (const auto& [input, reference] : bus_data.inputs)
          bus.addCurrent(static_cast<size_t>(input), source(reference));
      }

      for (const auto& model_data : data.inner_current_control)
      {
        auto& model  = component<Controller::InnerCurrentControl<ScalarT, IdxT>>(model_data.id);
        using Inputs = Controller::InnerCurrentControlInputs;
        std::array<SignalT*, static_cast<size_t>(Inputs::SIZE)> inputs{};
        for (size_t n = 0; n < inputs.size(); ++n)
          inputs[n] = &source(model_data.inputs.at(static_cast<Inputs>(n)));
        model.attachInput(inputs);
        for (const auto& [output, reference] : model_data.outputs)
          model.assignOutput(output, &signal(reference));
      }

      for (const auto& model_data : data.outer_voltage_control)
      {
        auto& model  = component<Controller::OuterVoltageControl<ScalarT, IdxT>>(model_data.id);
        using Inputs = Controller::OuterVoltageControlInputs;
        std::array<SignalT*, static_cast<size_t>(Inputs::SIZE)> inputs{};
        for (size_t n = 0; n < inputs.size(); ++n)
          inputs[n] = &source(model_data.inputs.at(static_cast<Inputs>(n)));
        model.attachInput(inputs);
        for (const auto& [output, reference] : model_data.outputs)
          model.assignOutput(output, &signal(reference));
      }
      for (const auto& model_data : data.outer_power_control)
      {
        auto& model  = component<Controller::OuterPowerControl<ScalarT, IdxT>>(model_data.id);
        using Inputs = Controller::OuterPowerControlInputs;
        std::array<SignalT*, static_cast<size_t>(Inputs::SIZE)> inputs{};
        for (size_t n = 0; n < inputs.size(); ++n)
          inputs[n] = &source(model_data.inputs.at(static_cast<Inputs>(n)));
        model.attachInput(inputs);
        for (const auto& [output, reference] : model_data.outputs)
          model.assignOutput(output, &signal(reference));
      }

      for (const auto& model_data : data.park)
      {
        auto& model = component<Park<ScalarT, IdxT>>(model_data.id);
        model.attachInput({&source(model_data.inputs.at(ParkInputs::u1)),
                           &source(model_data.inputs.at(ParkInputs::u2)),
                           &source(model_data.inputs.at(ParkInputs::u3))},
                          &source(model_data.inputs.at(ParkInputs::theta)));
        for (const auto& [output, reference] : model_data.outputs)
          model.assignOutput(output, &signal(reference));
      }

      for (const auto& model_data : data.pll)
      {
        auto& model = component<Pll<ScalarT, IdxT>>(model_data.id);
        model.attachInput({&source(model_data.inputs.at(PllInputs::va)),
                           &source(model_data.inputs.at(PllInputs::vb)),
                           &source(model_data.inputs.at(PllInputs::vc))});
        for (const auto& [output, id] : model_data.outputs)
          model.assignOutput(output, &signal(id));
      }

      for (const auto& model_data : data.pwm)
      {
        auto& model  = component<Controller::Pwm<ScalarT, IdxT>>(model_data.id);
        using Inputs = Controller::PwmInputs;
        if (!model_data.inputs.empty())
        {
          for (const auto key : {Inputs::ud, Inputs::uq, Inputs::vdc, Inputs::theta})
            if (!model_data.inputs.contains(key))
              throw std::invalid_argument("PWM \"" + model_data.id + "\" requires u, vdc, and theta inputs together");
          model.attachInput({&source(model_data.inputs.at(Inputs::ud)), &source(model_data.inputs.at(Inputs::uq))},
                            &source(model_data.inputs.at(Inputs::vdc)),
                            &source(model_data.inputs.at(Inputs::theta)));
        }
        for (const auto& [output, reference] : model_data.outputs)
          model.assignOutput(output, &signal(reference));
      }

      for (const auto& model_data : data.dc_link)
      {
        auto& model = component<Controller::DcLink<ScalarT, IdxT>>(model_data.id);
        model.attachInput(&source(model_data.inputs.at(Controller::DcLinkInputs::isrc)), &source(model_data.inputs.at(Controller::DcLinkInputs::idc)));
        for (const auto& [output, reference] : model_data.outputs)
          model.assignOutput(output, &signal(reference));
      }

      for (const auto& model_data : data.converter)
      {
        auto& model = component<Converter<ScalarT, IdxT>>(model_data.id);
        model.attachInput({&source(model_data.inputs.at(ConverterInputs::sa)),
                           &source(model_data.inputs.at(ConverterInputs::sb)),
                           &source(model_data.inputs.at(ConverterInputs::sc))},
                          &source(model_data.inputs.at(ConverterInputs::vdc)),
                          {&source(model_data.inputs.at(ConverterInputs::ia)),
                           &source(model_data.inputs.at(ConverterInputs::ib)),
                           &source(model_data.inputs.at(ConverterInputs::ic))});
        for (const auto& [output, reference] : model_data.outputs)
          model.assignOutput(output, &signal(reference));
      }

      for (const auto& filter_data : data.filter)
      {
        auto& filter = component<Filter<ScalarT, IdxT>>(filter_data.id);
        for (const auto& [output, reference] : filter_data.outputs)
          filter.assignOutput(output, &signal(reference));
        auto [bus, phases, voltage] = terminal(filter_data.inputs, FilterInputs::va);
        filter.attachInput(voltage, {&source(filter_data.inputs.at(FilterInputs::ea)), &source(filter_data.inputs.at(FilterInputs::eb)), &source(filter_data.inputs.at(FilterInputs::ec))});
        for (size_t p = 0; p < 3; ++p)
          bus->addCurrent(phases[p], filter.currentSignal(p));
      }

      for (const auto& source_data : data.voltage_source)
      {
        auto& source_model = component<VoltageSource<ScalarT, IdxT>>(source_data.id);
        for (const auto& [output, reference] : source_data.outputs)
          source_model.assignOutput(output, &signal(reference));
        auto [bus, phases, voltage] = terminal(source_data.inputs, VoltageSourceInputs::va);
        for (size_t p = 0; p < 3; ++p)
        {
          source_model.getSignals().attachSignal(static_cast<VoltageSourceExternalVariables>(p), voltage[p]);
          bus->addCurrent(phases[p], source_model.currentSignal(p));
        }
      }

      for (const auto& source_data : data.dependent_voltage_source)
      {
        auto& source_model = component<DependentVoltageSource<ScalarT, IdxT>>(source_data.id);
        for (const auto& [output, reference] : source_data.outputs)
          source_model.assignOutput(output, &signal(reference));
        auto [bus, phases, voltage] = terminal(source_data.inputs, DependentVoltageSourceInputs::va);
        for (size_t p = 0; p < 3; ++p)
        {
          source_model.getSignals().attachSignal(static_cast<DependentVoltageSourceExternalVariables>(p), voltage[p]);
          bus->addCurrent(phases[p], source_model.currentSignal(p));
        }

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
        auto& machine_model         = component<Machine<ScalarT, IdxT>>(machine_data.id);
        auto [bus, phases, voltage] = terminal(machine_data.inputs, MachineInputs::va);
        for (size_t p = 0; p < 3; ++p)
        {
          machine_model.getSignals().attachSignal(static_cast<MachineExternalVariables>(p), voltage[p]);
          bus->addCurrent(phases[p], machine_model.currentSignal(p));
        }

        for (const auto& [output, reference] : machine_data.outputs)
          machine_model.assignOutput(output, &signal(reference));
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

      for (const auto& source_data : data.regfma)
      {
        auto& source_model = component<Regfma<ScalarT, IdxT>>(source_data.id);
        for (const auto& [output, reference] : source_data.outputs)
          source_model.assignOutput(output, &signal(reference));
        auto [bus, phases, voltage] = terminal(source_data.inputs, RegfmaInputs::va);
        for (size_t p = 0; p < 3; ++p)
        {
          source_model.getSignals().attachSignal(static_cast<RegfmaExternalVariables>(p), voltage[p]);
          bus->addCurrent(phases[p], source_model.currentSignal(p));
        }
        for (const auto input : {RegfmaInputs::pref, RegfmaInputs::qref, RegfmaInputs::vref})
        {
          if (source_data.inputs.contains(input))
            source_model.getSignals().attachSignal(static_cast<RegfmaExternalVariables>(input),
                                                   &source(source_data.inputs.at(input)));
        }
      }

      for (const auto& line_data : data.line_lumped)
      {
        auto& line_model = component<LineLumped<ScalarT, IdxT>>(line_data.id);
        for (const auto& [output, reference] : line_data.outputs)
          line_model.assignOutput(output, &signal(reference));
        auto [bus1, phases1, voltage1] = terminal(line_data.inputs, LineLumpedInputs::v1a);
        auto [bus2, phases2, voltage2] = terminal(line_data.inputs, LineLumpedInputs::v2a);
        typename BusT::YDataT Y;
        using Parameter = LineLumpedParameters;
        if (line_data.Yp)
          Y = *line_data.Yp;
        else
        {
          if (line_data.parameters.contains(Parameter::Gp))
            Y.D = parameter<ABCMatrix<RealT>>(line_data, Parameter::Gp);
          if (line_data.parameters.contains(Parameter::Cp))
            Y.E = parameter<ABCMatrix<RealT>>(line_data, Parameter::Cp);
        }
        const RealT dx = parameter<RealT>(line_data, Parameter::dx, RealT{0});
        std::string name;
        for (char c : qualify(line_data.id))
          name += c == '/' ? "//" : c == '.' ? "/"
                                             : std::string(1, c);
        typename BusT::PhaseSignals i12, i21;
        for (size_t p = 0; p < 3; ++p)
        {
          i12[p] = &line_model.outputSignal(static_cast<LineLumpedOutputs>(p));
          i21[p] = &line_model.outputSignal(static_cast<LineLumpedOutputs>(3 + p));
        }
        for (size_t p = 0; p < 3; ++p)
        {
          bus1->addCurrent(phases1[p], *i21[p]);
          bus2->addCurrent(phases2[p], *i12[p]);
        }
        bus1->addShunt(name + "_1", Y, HALF<RealT> * dx, phases1);
        bus2->addShunt(name + "_2", Y, HALF<RealT> * dx, phases2);
        line_model.attachTerminal(0, voltage1);
        line_model.attachTerminal(1, voltage2);
      }

      for (const auto& line_data : data.line_distributed)
      {
        auto& line_model = component<LineDistributed<ScalarT, IdxT>>(line_data.id);
        for (const auto& [output, reference] : line_data.outputs)
          line_model.assignOutput(output, &signal(reference));
        auto [bus1, phases1, voltage1] = terminal(line_data.inputs, LineDistributedInputs::v1a);
        auto [bus2, phases2, voltage2] = terminal(line_data.inputs, LineDistributedInputs::v2a);
        std::string name;
        for (char c : qualify(line_data.id))
          name += c == '/' ? "//" : c == '.' ? "/"
                                             : std::string(1, c);
        std::array<BusT*, 2>                     buses{bus1, bus2};
        std::array<typename BusT::PhaseOrder, 2> phases{phases1, phases2};
        for (size_t e = 0; e < 2; ++e)
        {
          typename BusT::PhaseSignals incident, characteristic;
          for (size_t p = 0; p < 3; ++p)
            incident[p] = &line_model.incidentSignal(e, p);
          auto& norton = buses[e]->addNorton(name + "_" + std::to_string(e + 1), line_data.Yc, incident, ONE<RealT>, phases[e]);
          for (size_t p = 0; p < 3; ++p)
            characteristic[p] = &norton.outputSignal(p);
          line_model.attachTerminal(e, characteristic);
        }
      }

      for (const auto& load_data : data.loadz)
      {
        auto& load_model = component<LoadZ<ScalarT, IdxT>>(load_data.id);
        for (const auto& [output, reference] : load_data.outputs)
          load_model.getSignals().assignSignal(static_cast<LoadZInternalVariables>(output), &signal(reference));
        auto [bus, phases, voltage] = terminal(load_data.inputs, LoadZInputs::va);
        for (size_t p = 0; p < 3; ++p)
        {
          load_model.getSignals().attachSignal(static_cast<LoadZExternalVariables>(p), voltage[p]);
          bus->addCurrent(phases[p], load_model.currentSignal(p));
        }
      }

      for (const auto& stabilizer_data : data.ieeest)
      {
        auto& stabilizer_model = component<Controller::Ieeest<ScalarT, IdxT>>(stabilizer_data.id);
        if (stabilizer_data.inputs.contains(Controller::IeeestInputs::input))
          stabilizer_model.getSignals().template attachSignal<Controller::IeeestExternalVariables::U>(
              &source(stabilizer_data.inputs.at(Controller::IeeestInputs::input)));
        if (stabilizer_data.inputs.contains(Controller::IeeestInputs::speed))
          stabilizer_model.getSignals().template attachSignal<Controller::IeeestExternalVariables::OMEGA>(
              &source(stabilizer_data.inputs.at(Controller::IeeestInputs::speed)));
        if (stabilizer_data.inputs.contains(Controller::IeeestInputs::vct))
          stabilizer_model.getSignals().template attachSignal<Controller::IeeestExternalVariables::VCT>(
              &source(stabilizer_data.inputs.at(Controller::IeeestInputs::vct)));
        if (stabilizer_data.outputs.contains(Controller::IeeestOutputs::output))
          stabilizer_model.getSignals().template assignSignal<Controller::IeeestInternalVariables::VSS>(
              &signal(stabilizer_data.outputs.at(Controller::IeeestOutputs::output)));
      }

      for (const auto& governor_data : data.gastpti)
      {
        auto& governor_model = component<Controller::GastPti<ScalarT, IdxT>>(governor_data.id);
        if (governor_data.inputs.contains(Controller::GastPtiInputs::speed))
        {
          governor_model.getSignals().template attachSignal<Controller::GastPtiExternalVariables::OMEGA>(
              &source(governor_data.inputs.at(Controller::GastPtiInputs::speed)));
        }
        if (governor_data.inputs.contains(Controller::GastPtiInputs::pref))
        {
          governor_model.getSignals().template attachSignal<Controller::GastPtiExternalVariables::PREF>(
              &source(governor_data.inputs.at(Controller::GastPtiInputs::pref)));
        }
        if (governor_data.outputs.contains(Controller::GastPtiOutputs::pmech))
        {
          governor_model.getSignals().template assignSignal<Controller::GastPtiInternalVariables::PMECH>(
              &signal(governor_data.outputs.at(Controller::GastPtiOutputs::pmech)));
        }
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

      for (const auto& exciter_data : data.sexs_pti)
      {
        using Inputs        = Controller::SexsPtiInputs;
        using External      = Controller::SexsPtiExternalVariables;
        using Internal      = Controller::SexsPtiInternalVariables;
        auto& exciter_model = component<Controller::SexsPti<ScalarT, IdxT>>(exciter_data.id);
        auto& signals       = exciter_model.getSignals();
        signals.template attachSignal<External::VA>(&source(exciter_data.inputs.at(Inputs::va)));
        signals.template attachSignal<External::VB>(&source(exciter_data.inputs.at(Inputs::vb)));
        signals.template attachSignal<External::VC>(&source(exciter_data.inputs.at(Inputs::vc)));
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
        if (exciter_data.outputs.contains(Controller::SexsPtiOutputs::efd))
        {
          signals.template assignSignal<Internal::EFD>(
              &signal(exciter_data.outputs.at(Controller::SexsPtiOutputs::efd)));
        }
      }

      for (const auto& exciter_data : data.exciter)
      {
        using Inputs        = Controller::Ieeet1Inputs;
        using External      = Controller::Ieeet1ExternalVariables;
        using Internal      = Controller::Ieeet1InternalVariables;
        auto& exciter_model = component<Controller::Ieeet1<ScalarT, IdxT>>(exciter_data.id);
        auto& signals       = exciter_model.getSignals();
        signals.template attachSignal<External::VA>(&source(exciter_data.inputs.at(Inputs::va)));
        signals.template attachSignal<External::VB>(&source(exciter_data.inputs.at(Inputs::vb)));
        signals.template attachSignal<External::VC>(&source(exciter_data.inputs.at(Inputs::vc)));
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
        for (const auto& [output, reference] : switch_data.outputs)
          switch_model.getSignals().assignSignal(static_cast<SwitchInternalVariables>(output), &signal(reference));
        for (size_t end = 0; end < 2; ++end)
        {
          auto [bus, phases, voltage] = terminal(switch_data.inputs, static_cast<SwitchInputs>(3 * end));
          for (size_t p = 0; p < 3; ++p)
          {
            switch_model.getSignals().attachSignal(static_cast<SwitchExternalVariables>(3 * end + p), voltage[p]);
            bus->addCurrent(phases[p], switch_model.currentSignal(p), end == 0 ? -ONE<RealT> : ONE<RealT>);
          }
        }
      }

      for (const auto& transformer_data : data.transformer)
      {
        auto& transformer_model = component<Transformer<ScalarT, IdxT>>(transformer_data.id);
        for (const auto& [output, reference] : transformer_data.outputs)
          transformer_model.assignOutput(output, &signal(reference));
        for (size_t end = 0; end < 2; ++end)
        {
          auto [bus, phases, voltage] = terminal(transformer_data.inputs, static_cast<TransformerInputs>(3 * end));
          transformer_model.attachTerminal(end, voltage);
          for (size_t p = 0; p < 3; ++p)
            bus->addCurrent(phases[p], transformer_model.currentSignal(end, p));
        }
      }
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::SignalT*
    Container<scalar_type, index_type>::resolveOutput(std::string_view reference)
    {
      if (const auto found = signals_by_id_.find(reference); found != signals_by_id_.end())
        return found->second;
      const auto dot = reference.find('.');
      if (dot != std::string_view::npos)
      {
        auto&      child = component(reference.substr(0, dot));
        const auto name  = reference.substr(dot + 1);
        if (auto* container = dynamic_cast<Container*>(&child))
          return &container->outputSignal(name);
        if (auto* line = dynamic_cast<LineDistributed<ScalarT, IdxT>*>(&child))
        {
          const auto output = magic_enum::enum_cast<LineDistributedOutputs>(name);
          if (output && *output != LineDistributedOutputs::SIZE)
            return &line->outputSignal(*output);
        }
        if (auto* line = dynamic_cast<LineLumped<ScalarT, IdxT>*>(&child))
        {
          const auto output = magic_enum::enum_cast<LineLumpedOutputs>(name);
          if (output && *output != LineLumpedOutputs::SIZE)
            return &line->outputSignal(*output);
        }
      }
      throw std::invalid_argument("Unknown scalar signal: " + std::string(reference));
    }

  } // namespace EMT
} // namespace GridKit

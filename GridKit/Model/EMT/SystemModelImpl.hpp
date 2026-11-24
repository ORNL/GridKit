/**
 * @file SystemModelImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the EMT system model.
 *
 * The local variables, connections and indices of a component are defined in
 * `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <exception>
#include <format>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/ComponentLibrary.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct the buses, then the devices, and the monitor sinks of `data`.
     *
     * @throws std::invalid_argument for an unknown class, bus, device or name, or a duplicate id.
     *
     * @post Every device terminal is connected to its bus and every device
     *       input to what it reads, and `blocks_` lists each component before
     *       its submodels.
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data)
      : monitor_(std::make_unique<MonitorT>(this->time()))
    {
      this->id_               = "SystemModel";
      this->freq_system_base_ = data.freq_base;
      this->va_system_base_   = data.va_base;

      for (const auto& bus : data.buses)
      {
        if (!buses_.emplace(bus.id, &add(bus)).second)
        {
          throw std::invalid_argument(std::format("EMT::SystemModel: duplicate bus {}", bus.id));
        }
      }
      for (const auto& device : data.devices)
      {
        if (!devices_.emplace(device.id, &add(device)).second)
        {
          throw std::invalid_argument(std::format("EMT::SystemModel: duplicate device {}", device.id));
        }
      }

      // An input may read a device that the case lists later
      for (const auto& device : data.devices)
      {
        connect(device);
      }
      addSignals();

      for (const auto& sink : data.monitor_sink)
      {
        monitor_->addSink(sink);
      }
    }

    /**
     * @brief Add the library model of the class that `device` names.
     *
     * @param[in] device Case entry of a bus or a device.
     *
     * @throws std::invalid_argument if no library model has that class.
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::ComponentT&
    SystemModel<scalar_type, index_type>::add(const DeviceDataT& device)
    {
      ComponentT* added = nullptr;
      ComponentLibrary<ScalarT, IdxT>::forEach(
          [&]<typename ModelT>(std::type_identity<ModelT>)
          {
            if (ModelT::CLASS_NAME == device.device_class)
            {
              added = &addDevice<ModelT>(device);
            }
          });
      if (added == nullptr)
      {
        throw std::invalid_argument(std::format("{}: unknown class '{}'", device.id, device.device_class));
      }
      return *added;
    }

    /**
     * @brief Construct a `ModelT` from `device`, with its monitors.
     *
     * @tparam ModelT Library model of the entry.
     * @param[in] device Case entry of a bus or a device.
     *
     * @throws std::invalid_argument if the parameters are invalid or a
     *         monitor is unknown.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    ModelT& SystemModel<scalar_type, index_type>::addDevice(const DeviceDataT& device)
    {
      std::unique_ptr<ModelT> model;
      try
      {
        model = std::make_unique<ModelT>(device.params.template get<typename ModelT::ModelDataT>());
      }
      catch (const std::exception& error)
      {
        throw std::invalid_argument(std::format("{}: {}: {}", device.id, device.device_class, error.what()));
      }
      ComponentT& component = *model;
      component.id_         = device.id;

      ModelT& added = *model;
      addBlocks(component);
      for (const auto& name : device.mon)
      {
        addMonitor(ModelT::CLASS_NAME, added, name);
      }
      components_.push_back(std::move(model));
      return added;
    }

    /**
     * @brief Add the columns of monitor `name` of `model`.
     *
     * The name is an internal variable or an output of the model, or
     * `symbol.name` for one of its submodel `symbol`.
     *
     * @tparam ModelT Model of the device or of one of its submodels.
     * @param[in] device_class Class of the device.
     * @param[in] model Device, or submodel, that owns the monitor.
     * @param[in] name Monitor as the case names it.
     *
     * @throws std::invalid_argument if the name is unknown.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    void SystemModel<scalar_type, index_type>::addMonitor(std::string_view device_class, ModelT& model, const std::string& name)
    {
      if constexpr (requires { typename ModelT::SubmodelsT; })
      {
        const auto dot = name.find('.');
        if (dot != std::string::npos)
        {
          const auto symbol = Utilities::enumKey<typename ModelT::SubmodelsT>(name.substr(0, dot), "submodel", model.id());
          addMonitor(device_class, model.submodel(symbol), name.substr(dot + 1));
          return;
        }
      }

      if constexpr (!std::same_as<typename ModelT::OutputsT, Utilities::EmptyEnum>)
      {
        const auto output = Utilities::enum_cast<typename ModelT::OutputsT>(name);
        if (output.has_value())
        {
          monitored_.add(device_class, model, output.value());
          return;
        }
      }
      monitored_.add(device_class, model, Utilities::enumKey<typename ModelT::InternalVariablesT>(name, "monitor", model.id()));
    }

    /**
     * @brief Connect the terminals, inputs and outputs of `device` to what the case names.
     *
     * A terminal needs a bus. An input that the case leaves out is held,
     * unless the model has connected it, and an output that it leaves out
     * stays unconnected.
     *
     * @param[in] device Case entry of a device that `add` has constructed.
     *
     * @throws std::invalid_argument if a name, bus or device is unknown, a
     *         terminal is left out, or a connection does not fit what it names.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::connect(const DeviceDataT& device)
    {
      ComponentLibrary<ScalarT, IdxT>::visit(
          this->device(device.id),
          [&](auto& model)
          {
            using ModelT = std::remove_reference_t<decltype(model)>;

            // Connects each enumerator of `EnumT`, a `kind` of the model, to what the case names, and calls `absent` for one it leaves out
            const auto wire = [&]<typename EnumT>(std::type_identity<EnumT>,
                                                  std::string_view                          kind,
                                                  std::vector<Connection>&                  connections,
                                                  const std::map<std::string, std::string>& entries,
                                                  bool                                      reads_outputs,
                                                  auto                                      absent)
            {
              std::map<EnumT, std::string> references;
              for (const auto& [name, reference] : entries)
              {
                references[Utilities::enumKey<EnumT>(name, kind, device.id)] = reference;
              }
              for (const EnumT symbol : Utilities::enum_values<EnumT>())
              {
                Connection& connection = connections[static_cast<IdxT>(symbol)];
                const auto  reference  = references.find(symbol);
                if (reference == references.end())
                {
                  absent(connection, Utilities::enum_name(symbol));
                  continue;
                }
                if (connect(connection, reference->second, reads_outputs) != 0)
                {
                  throw std::invalid_argument(
                      std::format("{}: {} '{}' does not fit '{}'", device.id, kind, Utilities::enum_name(symbol), reference->second));
                }
              }
            };

            ComponentT& component = model;
            wire(std::type_identity<typename ModelT::TerminalsT>{}, "terminal", component.terminals_, device.terminals, false, [&device](Connection&, std::string_view name)
                 { throw std::invalid_argument(std::format("{}: terminal '{}' is not connected", device.id, name)); });
            wire(std::type_identity<typename ModelT::InputsT>{}, "input", component.inputs_, device.inputs, true, [this](Connection& input, std::string_view)
                 {
                   if (input.component == nullptr)
                   {
                     held_inputs_.push_back(&input);
                   } });
            wire(std::type_identity<typename ModelT::OutputsT>{}, "output", component.outputs_, device.outputs, false, [](Connection&, std::string_view) {});
          });
    }

    /**
     * @brief Connect `connection` to what `reference` names.
     *
     * A bus number names the voltages of the bus, and `owner.name` an
     * internal variable of device or bus `owner`: the connection is to those
     * variables. For an input, `owner.name` may instead name an output, which
     * the input reads through a signal that `addSignals` adds.
     *
     * @param[in,out] connection Terminal, input or output of a device.
     * @param[in] reference What the case names.
     * @param[in] reads_outputs Whether the connection is an input, which may read an output.
     * @return 0 on success, 1 if the connection does not fit what it names.
     *
     * @throws std::invalid_argument if the owner or the name is unknown.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::connect(Connection& connection, const std::string& reference, bool reads_outputs)
    {
      // Connects to the `size` local variables of `component` from `first`
      const auto variables = [&connection](ComponentT& component, IdxT first, IdxT size)
      {
        if (size != connection.indices.size())
        {
          return 1;
        }
        connection.component = &component;
        std::iota(connection.indices.begin(), connection.indices.end(), first);
        return 0;
      };

      const auto dot = reference.rfind('.');
      if (dot == std::string::npos)
      {
        ComponentT& target = bus(reference);
        return variables(target, 0, target.size_);
      }

      const std::string name   = reference.substr(dot + 1);
      ComponentT&       target = owner(reference.substr(0, dot));

      int status = 0;
      ComponentLibrary<ScalarT, IdxT>::visit(
          target,
          [&](auto& model)
          {
            using ModelT = std::remove_reference_t<decltype(model)>;

            if constexpr (!std::same_as<typename ModelT::OutputsT, Utilities::EmptyEnum>)
            {
              const auto named = Utilities::enum_cast<typename ModelT::OutputsT>(name);
              if (reads_outputs && named.has_value())
              {
                Connection& output = target.outputs_[static_cast<IdxT>(named.value())];
                if (output.indices.size() != connection.indices.size())
                {
                  status = 1;
                  return;
                }
                auto& [read, inputs] = signal_inputs_[reference];
                read                 = &output;
                inputs.push_back(&connection);
                return;
              }
            }

            const auto variable = Utilities::enumKey<typename ModelT::InternalVariablesT>(name, "variable", target.id());
            status              = variables(target, target.local(variable), model.size(variable));
          });
      return status;
    }

    /**
     * @brief Device `id`, or bus `id` when no device has that id.
     *
     * @throws std::invalid_argument if neither exists.
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::ComponentT& SystemModel<scalar_type, index_type>::owner(const std::string& id) const
    {
      const auto device = devices_.find(id);
      if (device != devices_.end())
      {
        return *device->second;
      }
      const auto bus = buses_.find(id);
      if (bus != buses_.end())
      {
        return *bus->second;
      }
      throw std::invalid_argument(std::format("EMT::SystemModel: unknown device or bus {}", id));
    }

    /**
     * @brief Add one `Signal` of held signals for the inputs that the case
     *        leaves out, and one for the outputs that inputs read.
     *
     * @throws std::invalid_argument if an output that an input reads is
     *         already connected.
     *
     * @post Every input of a device is connected.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addSignals()
    {
      // Adds `size` signals, held or not, as component `id`
      const auto block = [this](std::string_view id, IdxT size, bool held)
      {
        DeviceDataT data;
        data.device_class   = Signal<ScalarT, IdxT>::CLASS_NAME;
        data.id             = id;
        data.params["N"]    = size;
        data.params["held"] = held;
        return &add(data);
      };

      // Connects `connection` to the signals of `component` from `first`
      const auto assign = [](Connection& connection, ComponentT* component, IdxT first)
      {
        connection.component = component;
        std::iota(connection.indices.begin(), connection.indices.end(), first);
      };

      IdxT held = 0;
      for (const Connection* input : held_inputs_)
      {
        held += input->indices.size();
      }
      if (held > 0)
      {
        held_      = block("held", held, true);
        IdxT first = 0;
        for (Connection* input : held_inputs_)
        {
          assign(*input, held_, first);
          first += input->indices.size();
        }
      }

      IdxT signals = 0;
      for (const auto& [reference, read] : signal_inputs_)
      {
        signals += read.first->indices.size();
      }
      if (signals > 0)
      {
        signals_   = block("signal", signals, false);
        IdxT first = 0;
        for (const auto& [reference, read] : signal_inputs_)
        {
          const auto& [output, inputs] = read;
          if (output->component != nullptr)
          {
            throw std::invalid_argument(std::format("{}: output is already connected", reference));
          }

          assign(*output, signals_, first);
          for (Connection* input : inputs)
          {
            assign(*input, signals_, first);
          }
          first += output->indices.size();
        }
      }
    }

    /**
     * @brief Append `component`, then the submodels it has, to `blocks_`.
     *
     * Parents precede the submodels that read their variables. The id of a
     * submodel is `id.symbol`.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addBlocks(ComponentT& component)
    {
      blocks_.push_back(&component);
      if (const auto count = component.sizeRoots(); count != 0)
      {
        root_components_.push_back(&component);
        component_root_count_ += count;
      }
      for (const auto& submodel : component.submodels_)
      {
        if (submodel)
        {
          submodel->id_ = std::format("{}.{}", component.id_, submodel->id_);
          addBlocks(*submodel);
        }
      }
    }

    /**
     * @brief Check every component and its connections within this system.
     *
     * A terminal and an input are connected, and an output may be. A
     * connection to an external variable of another component requires that
     * component earlier in `blocks_`.
     *
     * @return The failed checks, each naming its component.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks SystemModel<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks;

      std::map<const ComponentT*, IdxT> position;
      for (IdxT b = 0; b < blocks_.size(); ++b)
      {
        position.emplace(blocks_[b], b);
      }

      for (IdxT b = 0; b < blocks_.size(); ++b)
      {
        const ComponentT* component        = blocks_[b];
        const auto        component_checks = component->verify();
        for (const auto& error : component_checks.errors())
        {
          checks.fail(std::format("{}: {}", component->id(), error));
        }

        const auto connected = [&](const std::vector<Connection>& connections, std::string_view kind, bool required)
        {
          for (IdxT c = 0; c < connections.size(); ++c)
          {
            const Connection& connection = connections[c];
            if (connection.component == nullptr && !required)
            {
              continue;
            }

            const auto found = position.find(connection.component);
            if (found == position.end())
            {
              checks.fail(std::format("{}: {} {} is not connected to this system", component->id(), kind, c));
              continue;
            }

            // An external variable of another component resolves through that component first
            const bool forwarded = std::ranges::any_of(connection.indices, [&connection](IdxT i)
                                                       { return i >= connection.component->size_; });
            checks.check(!forwarded || found->second < b,
                         std::format("{}: {} {} is connected to a later component", component->id(), kind, c));
          }
        };
        connected(component->terminals_, "terminal", true);
        connected(component->inputs_, "input", true);
        connected(component->outputs_, "output", false);
      }

      return checks;
    }

    /**
     * @brief Verify the system, compose the system vectors, allocate the
     *        Jacobian and start the monitor.
     *
     * @return 0 on success, 1 if a check fails or an allocation fails.
     *
     * @post Each component views its block of the system vectors and has
     *       resolved its indices.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocate()
    {
      const auto checks = verify();
      for (const auto& error : checks.errors())
      {
        Log::error() << error << '\n';
      }
      if (!checks.passed())
      {
        return 1;
      }

      if (compose() != 0 || allocateDerivatives() != 0)
      {
        Log::error() << "EMT::SystemModel: failed to allocate the system vectors and derivatives\n";
        return 1;
      }

      monitor_->addMonitor(&monitored_);
      monitor_->start();
      return 0;
    }

    /**
     * @brief Mark the values and derivatives that the state has set.
     *
     * Components initialize their discrete modes from the supplied state.
     *
     * @return 0.
     *
     * @pre `allocate` has succeeded.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      for (ComponentT* component : blocks_)
      {
        if (const int status = component->initialize(); status != 0)
        {
          return status;
        }
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return 0;
    }

    /**
     * @brief Set the absolute tolerance of every variable, in block order.
     *
     * Buses come first and parents before submodels, so each model reads the
     * tolerances of its external values. A signal takes the tolerance that
     * the model of its output sets, which defaults to `rel_tol`.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0 on success, or the status of the first component that fails.
     *
     * @pre `allocate` has succeeded.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      for (ComponentT* component : blocks_)
      {
        std::ranges::fill(component->output_tol_, static_cast<ScalarT>(rel_tol));
        const int status = component->setAbsoluteTolerance(rel_tol);
        if (status != 0)
        {
          return status;
        }
      }

      // A signal takes the tolerance of the output it follows
      for (const ComponentT* component : blocks_)
      {
        IdxT entry = 0;
        for (const Connection& output : component->outputs_)
        {
          if (signals_ != nullptr && output.component == signals_)
          {
            ScalarT* tolerance = signals_->abs_tol_.getData();
            for (IdxT k = 0; k < output.indices.size(); ++k)
            {
              tolerance[output.indices[k]] = component->output_tol_[entry + k];
            }
          }
          entry += output.indices.size();
        }
      }
      abs_tol_.setDataUpdated();
      return 0;
    }

    /**
     * @brief Zero the system residual, then add the residual of every component.
     *
     * @return 0 on success, or the status of the first component that fails.
     *
     * @pre `allocate` has succeeded.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateResidual()
    {
      f_.setToZero();
      for (ComponentT* component : blocks_)
      {
        const int status = component->evaluateResidual();
        if (status != 0)
        {
          return status;
        }
      }
      f_.setDataUpdated();
      return 0;
    }

    /**
     * @brief Number of system variables, zero before `allocate`.
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::IdxT SystemModel<scalar_type, index_type>::size()
    {
      return y_.getSize();
    }

    /**
     * @brief Number of nonzeros of the system Jacobian, zero before `allocate`.
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::IdxT SystemModel<scalar_type, index_type>::nnz()
    {
      IdxT nnz = 0;
      if (jacobian_)
      {
        nnz = jacobian_->getNnz();
      }
      return nnz;
    }

    /**
     * @brief System Jacobian in CSR format, null before `allocate`.
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::CsrMatrixT* SystemModel<scalar_type, index_type>::getCsrJacobian() const
    {
      return jacobian_.get();
    }

    /**
     * @brief True if the monitor has entries to write.
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::monitoring() const
    {
      return !monitor_->empty();
    }

    /**
     * @brief Write the monitored variables to the monitor sinks.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::printMonitoredVariables() const
    {
      monitored_.evaluateOutputs();
      monitor_->print();
    }

    /**
     * @brief Start the monitor output.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::startMonitor()
    {
      monitor_->start();
    }

    /**
     * @brief Stop the monitor output.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::stopMonitor()
    {
      monitor_->stop();
    }

    /**
     * @brief Bus of case number `number`.
     *
     * @throws std::invalid_argument if the system has no such bus.
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::ComponentT&
    SystemModel<scalar_type, index_type>::bus(const std::string& number) const
    {
      const auto found = buses_.find(number);
      if (found == buses_.end())
      {
        throw std::invalid_argument(std::format("EMT::SystemModel: unknown bus {}", number));
      }
      return *found->second;
    }

    /**
     * @brief Device of case id `id`.
     *
     * @throws std::invalid_argument if the system has no such device.
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::ComponentT&
    SystemModel<scalar_type, index_type>::device(const std::string& id) const
    {
      const auto found = devices_.find(id);
      if (found == devices_.end())
      {
        throw std::invalid_argument(std::format("EMT::SystemModel: unknown device {}", id));
      }
      return *found->second;
    }

    /**
     * @brief Size the system vectors, bind each component to its block, then
     *        resolve the external variables.
     *
     * @return 0 on success, 1 if a vector cannot be sized or bound.
     *
     * @pre `verify` has passed, so each component follows the components
     *      whose external variables it is connected to.
     * @post The system values and derivatives are zero, and every component
     *       has resolved its indices.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::compose()
    {
      IdxT size = 0;
      for (const ComponentT* component : blocks_)
      {
        size += component->size_;
      }

      if (y_.resize(size) != 0
          || yp_.resize(size) != 0
          || f_.resize(size) != 0
          || abs_tol_.resize(size) != 0
          || y_.setToZero() != 0
          || yp_.setToZero() != 0)
      {
        return 1;
      }
      tag_.assign(size, false);

      IdxT offset = 0;
      for (ComponentT* component : blocks_)
      {
        if (bind(*component, offset) != 0)
        {
          return 1;
        }
        offset += component->size_;
      }

      // An external variable of a parent resolves through the parent, which comes first
      for (ComponentT* component : blocks_)
      {
        component->setIndices();
      }
      return 0;
    }

    /**
     * @brief Bind `component` to its block of the system vectors.
     *
     * @param[in] component Component that owns the block.
     * @param[in] offset System index of its first internal variable.
     * @return 0 on success, 1 if a view cannot be set.
     *
     * @pre The system vectors are sized.
     * @post `component` views its block, holds the system index of each
     *       internal variable and residual, shares the system time, and holds
     *       the system bases.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::bind(ComponentT& component, IdxT offset)
    {
      const auto view = [](auto& local, auto& global, IdxT start, IdxT count)
      { return local.setData(global.getData() + start, count) != 0; };

      const IdxT size = component.size_;
      if (view(component.y_, y_, offset, size)
          || view(component.yp_, yp_, offset, size)
          || view(component.f_, f_, offset, size)
          || view(component.abs_tol_, abs_tol_, offset, size))
      {
        return 1;
      }

      // A connection to an internal variable of a later component resolves once every component is bound
      IdxT* indices = component.variable_indices_.getData();
      std::iota(indices, indices + size, offset);
      std::iota(component.residual_indices_.data(), component.residual_indices_.data() + size, offset);

      component.system_y_       = y_.getData();
      component.system_yp_      = yp_.getData();
      component.system_f_       = f_.getData();
      component.system_abs_tol_ = abs_tol_.getData();

      component.clock_            = this->clock_;
      component.freq_system_base_ = this->freq_system_base_;
      component.va_system_base_   = this->va_system_base_;
      return 0;
    }
  } // namespace EMT
} // namespace GridKit

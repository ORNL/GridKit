/**
 * @file SystemModelImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the EMT system model.
 */

#pragma once

#include <algorithm>
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
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data)
      : monitor_(std::make_unique<MonitorT>(time_))
    {
      this->id_             = "SystemModel";
      this->va_system_base_ = data.va_base;

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

      for (const auto& sink : data.monitor_sink)
      {
        monitor_->addSink(sink);
      }
    }

    /// Library model of the entry's class
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
        throw std::invalid_argument(std::format("{} {}: unknown class", device.device_class, device.id));
      }
      return *added;
    }

    /// Model of `device`, connected at each port to the bus it names
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    ModelT& SystemModel<scalar_type, index_type>::addDevice(const DeviceDataT& device)
    {
      using PortsT = typename ModelT::PortsT;

      std::unique_ptr<ModelT> model;
      try
      {
        model = std::make_unique<ModelT>(device.params.template get<typename ModelT::ModelDataT>());
      }
      catch (const std::exception& error)
      {
        throw std::invalid_argument(std::format("{} {}: {}", device.device_class, device.id, error.what()));
      }
      ComponentT& component = *model;
      component.id_         = device.id;

      std::map<PortsT, IdxT> ports;
      for (const auto& [name, number] : device.ports)
      {
        ports[Utilities::enumKey<PortsT>(name, "port", device.id)] = number;
      }
      for (const PortsT port : Utilities::enum_values<PortsT>())
      {
        const auto number = ports.find(port);
        if (number == ports.end())
        {
          throw std::invalid_argument(std::format("{}: port '{}' has no bus", device.id, Utilities::enum_name(port)));
        }

        ComponentT& target = bus(std::to_string(number->second));
        if (component.ports_[static_cast<IdxT>(port)].indices.size() != target.size_ || component.connect(port, target) != 0)
        {
          throw std::invalid_argument(
              std::format("{}: port '{}' does not fit bus {}", device.id, Utilities::enum_name(port), number->second));
        }
      }

      for (const auto& name : device.mon)
      {
        monitored_.add(ModelT::CLASS_NAME, component, Utilities::enumKey<typename ModelT::InternalVariablesT>(name, "variable", device.id));
      }

      if constexpr (requires { typename ModelT::SubmodelsT; })
      {
        for (const auto symbol : Utilities::enum_values<typename ModelT::SubmodelsT>())
        {
          ComponentT& submodel = model->submodel(symbol);
          submodel.id_         = std::format("{}.{}", device.id, Utilities::enum_name(symbol));
        }
      }

      ModelT& added = *model;
      addBlocks(component);
      components_.push_back(std::move(model));
      return added;
    }

    /// Parents precede the submodels that read their variables
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addBlocks(ComponentT& component)
    {
      blocks_.push_back(&component);
      for (const auto& submodel : component.submodels_)
      {
        addBlocks(*submodel);
      }
    }

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

        for (IdxT p = 0; p < component->ports_.size(); ++p)
        {
          const Port& port  = component->ports_[p];
          const auto  found = position.find(port.component);
          if (found == position.end())
          {
            checks.fail(std::format("{}: port {} is not connected to this system", component->id(), p));
            continue;
          }

          // An external variable on another component's port resolves through that component first
          const bool forwarded = std::ranges::any_of(port.indices, [&port](IdxT i)
                                                     { return i >= port.component->size_; });
          checks.check(!forwarded || found->second < b,
                       std::format("{}: port {} reads a port of a later component", component->id(), p));
        }
      }

      return checks;
    }

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
        return 1;
      }

      monitor_->addMonitor(&monitored_);
      monitor_->start();
      return 0;
    }

    /// The initial values come from the state
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      for (ComponentT* component : blocks_)
      {
        const int status = component->initialize();
        if (status != 0)
        {
          Log::error() << component->id() << ": initialization failed\n";
          return status;
        }
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();
      return 0;
    }

    /// Buses first and parents before submodels, so each model reads the tolerances at its ports
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      for (ComponentT* component : blocks_)
      {
        const int status = component->setAbsoluteTolerance(rel_tol);
        if (status != 0)
        {
          return status;
        }
      }
      abs_tol_.setDataUpdated();
      return 0;
    }

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

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::updateTime(RealT t, RealT alpha)
    {
      time_        = t;
      this->alpha_ = alpha;
      for (ComponentT* component : blocks_)
      {
        component->updateTime(t, alpha);
      }
    }

    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::IdxT SystemModel<scalar_type, index_type>::size()
    {
      return y_.getSize();
    }

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

    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::CsrMatrixT* SystemModel<scalar_type, index_type>::getCsrJacobian() const
    {
      return jacobian_.get();
    }

    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::monitoring() const
    {
      return !monitor_->empty();
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::printMonitoredVariables() const
    {
      monitor_->print();
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::startMonitor()
    {
      monitor_->start();
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::stopMonitor()
    {
      monitor_->stop();
    }

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

    /// Size the system vectors, bind each component to its blocks, then resolve the external variables
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

      // An external variable on a parent's port resolves through the parent, which comes first
      for (ComponentT* component : blocks_)
      {
        component->setIndices();
      }
      return 0;
    }

    /// Views of the system vectors from `offset`
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

      IdxT* indices = component.variable_indices_.getData();
      std::iota(indices, indices + size, offset);

      component.system_y_       = y_.getData();
      component.system_yp_      = yp_.getData();
      component.system_f_       = f_.getData();
      component.system_abs_tol_ = abs_tol_.getData();

      component.va_system_base_ = this->va_system_base_;
      return 0;
    }
  } // namespace EMT
} // namespace GridKit

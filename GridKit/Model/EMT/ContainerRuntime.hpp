#pragma once

#include <algorithm>

#include <GridKit/Model/EMT/Container.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Container<scalar_type, index_type>::Container()
    {
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
                                                       SignalT*         value)
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
    typename Container<scalar_type, index_type>::SignalT*
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
      return *value;
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::inputSignal(std::string_view name) const
    {
      const auto value = inputEndpoint(name);
      return *value;
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
    typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::outputSignal(std::string_view name)
    {
      const auto found = outputs_.find(name);
      if (found == outputs_.end())
      {
        throw std::invalid_argument("Container output \"" + std::string(name)
                                    + "\" is not a scalar signal");
      }
      return *found->second;
    }

    template <typename scalar_type, typename index_type>
    const typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::outputSignal(std::string_view name) const
    {
      const auto found = outputs_.find(name);
      if (found == outputs_.end())
      {
        throw std::invalid_argument("Container output \"" + std::string(name)
                                    + "\" is not a scalar signal");
      }
      return *found->second;
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::SignalT&
    Container<scalar_type, index_type>::source(std::string_view reference)
    {
      const auto value = endpoint(reference);
      return *value;
    }

    template <typename scalar_type, typename index_type>
    typename Container<scalar_type, index_type>::SignalT*
    Container<scalar_type, index_type>::endpoint(std::string_view reference)
    {
      if (const auto found = inputs_.find(reference); found != inputs_.end())
        return found->second;
      return resolveOutput(reference);
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
        if (!value->hasProducer())
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
    int Container<scalar_type, index_type>::initialize(const std::map<std::string, std::map<std::string, RealT>>& state, RealT omega)
    {
      std::map<ComponentT*, std::string> paths;
      auto                               collect = [&](auto&& self, Container& scope, const std::string& prefix) -> void
      {
        for (const auto& [name, child] : scope.children_by_id_)
        {
          const auto path = prefix + name;
          paths.emplace(child, path);
          if (auto* container = dynamic_cast<Container*>(child))
          {
            self(self, *container, path + ".");
            if (auto* owner = container->initialStateComponent())
              paths[owner] = path;
          }
        }
      };
      collect(collect, *this, "");
      std::vector<ComponentT*> leaves;
      forEachComponent([&](ComponentT& component)
                       {
                         if (dynamic_cast<Container*>(&component) == nullptr)
                           leaves.push_back(&component); });
      std::set<std::string> initial_paths;
      for (auto* leaf : leaves)
        initial_paths.insert(paths.at(leaf));
      for (const auto& [path, values] : state)
        if (!initial_paths.contains(path))
          throw std::invalid_argument("Unknown initial state path: " + path);

      typename ComponentT::InitialStateT initial(omega);
      const std::map<std::string, RealT> empty;
      for (auto* leaf : leaves)
      {
        const auto& path  = paths.at(leaf);
        const auto  entry = state.find(path);
        initial.add(*leaf, path, entry == state.end() ? empty : entry->second, leaf->initializationPorts());
      }
      const int status = initial.initialize();
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return status;
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::setDifferentialTags(const std::set<size_t>& columns)
    {
      ComponentT::setDifferentialTags(columns);
      for (auto& child : children_)
        child->setDifferentialTags(columns);
    }

    template <typename scalar_type, typename index_type>
    std::string Container<scalar_type, index_type>::describeDaeIndex(IdxT index) const
    {
      for (const auto& [name, child] : children_by_id_)
      {
        const auto& indices = child->getVariableIndices();
        const auto  found   = std::find(indices.begin(), indices.end(), index);
        if (found == indices.end())
          continue;
        if (const auto* scope = dynamic_cast<const Container*>(child))
          return name + "." + scope->describeDaeIndex(index);
        return name + "[" + std::to_string(found - indices.begin()) + "]";
      }
      return "[" + std::to_string(index) + "]";
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
    int Container<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      size_t required = 0;
      for (auto& child : children_)
      {
        const int status = child->evaluateJacobian(y_scale, yp_scale);
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

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::resetHistory()
    {
      ComponentT::resetHistory();
      for (auto& child : children_)
        child->resetHistory();
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::acceptStep(RealT time)
    {
      ComponentT::acceptStep(time);
      for (auto& child : children_)
        child->acceptStep(time);
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::forEachChild(const std::function<void(const ComponentT&)>& visitor) const
    {
      for (const auto& child : children_)
        visitor(*child);
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::beginDiscontinuity(RealT time)
    {
      ComponentT::beginDiscontinuity(time);
      for (auto& child : children_)
        child->beginDiscontinuity(time);
    }

    template <typename scalar_type, typename index_type>
    auto Container<scalar_type, index_type>::maximumStepSize() const -> RealT
    {
      RealT step = ComponentT::maximumStepSize();
      for (const auto& child : children_)
      {
        const auto limit = child->maximumStepSize();
        if (!(limit > RealT{0}))
          throw std::invalid_argument("EMT history step bound must be positive");
        step = std::min(step, limit);
      }
      return step;
    }
  } // namespace EMT
} // namespace GridKit

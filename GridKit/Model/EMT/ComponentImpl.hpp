/**
 * @file ComponentImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the EMT component base class.
 */

#pragma once

#include <algorithm>
#include <new>
#include <numeric>
#include <utility>

#include <GridKit/Model/EMT/Component.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Component<scalar_type, index_type>::Component()
    {
      // Sized even when empty, since an unsized vector reports stale data
      if (jacobian_entries_.resize(0) != 0 || dfdy_.resize(0) != 0 || dfdyp_.resize(0) != 0)
      {
        Log::error() << "failed to allocate the derivative entries\n";
        throw std::bad_alloc{};
      }
    }

    template <typename scalar_type, typename index_type>
    template <typename PortT>
    int Component<scalar_type, index_type>::connect(PortT port, Component& component, IdxT offset)
    {
      const auto index = static_cast<IdxT>(port);
      if (index >= ports_.size())
      {
        return 1;
      }

      std::vector<IdxT> indices(ports_[index].indices.size());
      std::iota(indices.begin(), indices.end(), offset);
      return connect(port, component, std::move(indices));
    }

    template <typename scalar_type, typename index_type>
    template <typename PortT>
    int Component<scalar_type, index_type>::connect(PortT port, Component& component, std::vector<IdxT> indices)
    {
      const auto index = static_cast<IdxT>(port);
      if (index >= ports_.size() || indices.size() != ports_[index].indices.size())
      {
        return 1;
      }

      const IdxT local_size = component.variable_indices_.getSize();
      if (std::ranges::any_of(indices, [local_size](IdxT i)
                              { return i >= local_size; }))
      {
        return 1;
      }

      ports_[index] = {&component, std::move(indices)};
      return 0;
    }

    template <typename scalar_type, typename index_type>
    template <typename PortT>
    typename Component<scalar_type, index_type>::ScalarT
    Component<scalar_type, index_type>::absoluteTolerance(PortT port, IdxT k) const
    {
      return system_abs_tol_[variable_indices_.getData()[local(port, k)]];
    }

    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateJacobian()
    {
      const RealT* dfdy   = dfdy_.getData();
      const RealT* dfdyp  = dfdyp_.getData();
      RealT*       values = jacobian_entries_.values.getData();
      for (IdxT k = 0; k < jacobian_entries_.nnz; ++k)
      {
        values[k] = dfdy[k] + alpha_ * dfdyp[k];
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::setLayout(const std::vector<IdxT>& variables, const std::vector<IdxT>& ports)
    {
      offsets_.assign(1, 0);
      for (const IdxT variable : variables)
      {
        offsets_.push_back(offsets_.back() + variable);
      }
      size_ = offsets_.back();

      ports_.clear();
      IdxT local_size = size_;
      for (const IdxT port : ports)
      {
        ports_.push_back({nullptr, std::vector<IdxT>(port)});
        local_size += port;
      }

      if (variable_indices_.resize(local_size) != 0
          || local_y_.resize(local_size) != 0
          || local_yp_.resize(local_size) != 0
          || local_f_.resize(local_size) != 0)
      {
        Log::error() << id_ << ": failed to allocate the local variables\n";
        throw std::bad_alloc{};
      }
    }

    template <typename scalar_type, typename index_type>
    template <Utilities::SizedEnum PortT>
    typename Component<scalar_type, index_type>::IdxT Component<scalar_type, index_type>::local(PortT port, IdxT k) const
    {
      return size_ + external(port) + k;
    }

    template <typename scalar_type, typename index_type>
    template <Utilities::SizedEnum PortT>
    typename Component<scalar_type, index_type>::IdxT Component<scalar_type, index_type>::external(PortT port) const
    {
      IdxT index = 0;
      for (IdxT p = 0; p < static_cast<IdxT>(port); ++p)
      {
        index += ports_[p].indices.size();
      }
      return index;
    }

    template <typename scalar_type, typename index_type>
    template <typename ModelT, Utilities::SizedEnum SymbolT, typename... Args>
    ModelT& Component<scalar_type, index_type>::addSubmodel(SymbolT symbol, Args&&... args)
    {
      const auto index = static_cast<IdxT>(symbol);
      if (submodels_.size() <= index)
      {
        submodels_.resize(index + 1);
      }

      auto    submodel  = std::make_unique<ModelT>(std::forward<Args>(args)...);
      ModelT& added     = *submodel;
      submodels_[index] = std::move(submodel);
      return added;
    }

    template <typename scalar_type, typename index_type>
    template <typename ModelT, Utilities::SizedEnum SymbolT>
    ModelT& Component<scalar_type, index_type>::submodel(SymbolT symbol)
    {
      return static_cast<ModelT&>(*submodels_[static_cast<IdxT>(symbol)]);
    }

    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    int Component<scalar_type, index_type>::evaluateResidual(const ModelT& model)
    {
      ScalarT* y  = local_y_.getData();
      ScalarT* yp = local_yp_.getData();
      ScalarT* f  = local_f_.getData();

      gather(y, yp);
      local_f_.setToZero();
      model.internalResidual(y, yp, y + size_, yp + size_, f);
      model.externalResidual(y, yp, y + size_, yp + size_, f + size_);
      scatter(f);
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::gather(ScalarT* y, ScalarT* yp) const
    {
      const IdxT* indices = variable_indices_.getData();
      for (IdxT j = 0; j < variable_indices_.getSize(); ++j)
      {
        y[j]  = system_y_[indices[j]];
        yp[j] = system_yp_[indices[j]];
      }
    }

    /// Devices at one bus add to the same residuals, so a parallel scatter must synchronize
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::scatter(const ScalarT* f)
    {
      const IdxT* indices = variable_indices_.getData();
      for (IdxT j = 0; j < variable_indices_.getSize(); ++j)
      {
        system_f_[indices[j]] += f[j];
      }
    }

    /// An external variable on another component's port resolves through that component's indices
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::setIndices()
    {
      IdxT* indices = variable_indices_.getData() + size_;
      for (const Port& port : ports_)
      {
        for (const IdxT i : port.indices)
        {
          *indices = port.component->variable_indices_.getData()[i];
          ++indices;
        }
      }
    }
  } // namespace EMT
} // namespace GridKit

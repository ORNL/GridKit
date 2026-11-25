/**
 * @file ComponentImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the EMT component base class.
 *
 * Local values of a component, in order:
 * - internal variables, `size_` entries, variable `v` starting at `local(v)`
 * - external values: those of each terminal, then those of each input
 *
 * Local residuals of a component, in order:
 * - those of the internal variables
 * - external residuals: those of each terminal, then those of each output
 *
 * `variable_indices_` maps each local value to its system variable, and
 * `residual_indices_` each local residual to its system residual.
 *
 * Kernels of a model take `(y, yp, w, wp, f)`:
 * - `y`, `yp`: internal values and derivatives, indexed by `local(v)`
 * - `w`, `wp`: external values and derivatives, indexed by
 *   `terminalOffset(t) + k` or `inputOffset(i) + k`
 * - `f`: the block the kernel writes. `internalResidual` writes the internal
 *   residuals. `externalResidual`, of a model with terminals, writes what the
 *   model adds at its terminals, indexed by `terminalOffset(t) + k`.
 *
 * A model with outputs provides `output(o, y, yp, w, wp, g)`, which writes the
 * entries of output `o` to `g`, from its first. An output is evaluated with
 * the residual only when it is connected, and otherwise only by
 * `evaluateOutputs`, for a monitor.
 *
 * A model forwards `tagDifferentiable`, `evaluateOutputs` if it has outputs,
 * and `evaluateResidual` and `evaluateJacobian` if it is not linear, to this
 * class from the file that instantiates it. A linear model forwards neither,
 * since this class evaluates both from the partials that `tagDifferentiable`
 * stores. A model that is linear but for a term in time forwards
 * `evaluateResidual` only.
 *
 * The scalar type `T` of a kernel is `ScalarT` for the residual and a
 * dependency-tracking variable for the partials.
 *
 * A component adds to the system residual, so the system zeroes the residual
 * before evaluating its components.
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <format>
#include <new>
#include <numeric>
#include <stdexcept>
#include <utility>

#include <GridKit/Model/EMT/Component.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct an empty component.
     *
     * A model sizes itself with `setLayout`.
     *
     * @throws std::bad_alloc if the derivative entries cannot be allocated.
     */
    template <typename scalar_type, typename index_type>
    Component<scalar_type, index_type>::Component()
    {
      // Sized even when empty, since an unsized vector reports stale data
      if (jacobian_entries_.resize(0) != 0 || dfdy_.resize(0) != 0 || dfdyp_.resize(0) != 0)
      {
        Log::error() << "EMT::Component: failed to allocate the derivative entries\n";
        throw std::bad_alloc{};
      }
    }

    /**
     * @brief Connect `input` to the local values of `component` at `indices`.
     *
     * @tparam InputT Input enumeration of the model.
     * @param[in] input Input to connect.
     * @param[in] component Component whose values the input reads.
     * @param[in] indices Local index in `component` of each of those values.
     * @return 0 on success, 1 if `indices` has the wrong size or an index is out of range.
     *
     * @pre `setLayout` has sized both components.
     */
    template <typename scalar_type, typename index_type>
    template <Utilities::SizedEnum InputT>
    int Component<scalar_type, index_type>::connectInput(InputT input, Component& component, std::vector<IdxT> indices)
    {
      return connect(inputs_[static_cast<IdxT>(input)], component, std::move(indices), component.variable_indices_.getSize());
    }

    /**
     * @brief Connect `output` to the local residuals of `component` at `indices`.
     *
     * @tparam OutputT Output enumeration of the model.
     * @param[in] output Output to connect.
     * @param[in] component Component whose residuals the output is added to.
     * @param[in] indices Local index in `component` of each of those residuals.
     * @return 0 on success, 1 if `indices` has the wrong size or an index is out of range.
     *
     * @pre `setLayout` has sized both components.
     */
    template <typename scalar_type, typename index_type>
    template <Utilities::SizedEnum OutputT>
    int Component<scalar_type, index_type>::connectOutput(OutputT output, Component& component, std::vector<IdxT> indices)
    {
      return connect(outputs_[static_cast<IdxT>(output)], component, std::move(indices), component.residual_indices_.size());
    }

    /**
     * @brief Set `connection` to the local variables of `component` at `indices`.
     *
     * @param[in,out] connection Connection of a terminal, an input or an output.
     * @param[in] component Component to connect to.
     * @param[in] indices Local index in `component` of each variable.
     * @param[in] limit Number of local variables of `component` that the connection may name.
     * @return 0 on success, 1 if `indices` has the wrong size or an index reaches `limit`.
     */
    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::connect(Connection& connection, Component& component, std::vector<IdxT> indices, IdxT limit)
    {
      if (indices.size() != connection.indices.size())
      {
        return 1;
      }
      if (std::ranges::any_of(indices, [limit](IdxT i)
                              { return i >= limit; }))
      {
        return 1;
      }

      connection = {&component, std::move(indices)};
      return 0;
    }

    /**
     * @brief Absolute tolerance of external value `k`.
     *
     * @param[in] k Index among the external values, `terminalOffset(t) + k` or `inputOffset(i) + k`.
     * @return Tolerance set by the component that owns the variable.
     *
     * @pre The system has resolved the indices of this component, and the
     *      owner has set its tolerance.
     */
    template <typename scalar_type, typename index_type>
    typename Component<scalar_type, index_type>::ScalarT Component<scalar_type, index_type>::externalTolerance(IdxT k) const
    {
      return system_abs_tol_[external_indices_[k]];
    }

    /**
     * @brief Add the residual of a linear model to the system residual.
     *
     * The residual of a linear model is its constant partials dF/dy and
     * dF/dy' applied to the values and the derivatives, so no kernel is
     * evaluated. A model with any other term forwards `evaluateResidual`
     * to its kernels.
     *
     * @return 0.
     *
     * @pre `tagDifferentiable` has stored the partials.
     */
    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateResidual()
    {
      const IdxT*  rows  = jacobian_entries_.rows.getData();
      const IdxT*  cols  = jacobian_entries_.cols.getData();
      const RealT* dfdy  = dfdy_.getData();
      const RealT* dfdyp = dfdyp_.getData();
      for (IdxT k = 0; k < jacobian_entries_.nnz; ++k)
      {
        system_f_[rows[k]] += dfdy[k] * system_y_[cols[k]] + dfdyp[k] * system_yp_[cols[k]];
      }
      return 0;
    }

    /**
     * @brief Set the Jacobian entries to dF/dy + alpha dF/dy'.
     *
     * @return 0.
     *
     * @pre `tagDifferentiable` has stored the partials of the entries.
     */
    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateJacobian()
    {
      const RealT* dfdy   = dfdy_.getData();
      const RealT* dfdyp  = dfdyp_.getData();
      RealT*       values = jacobian_entries_.values.getData();
      for (IdxT k = 0; k < jacobian_entries_.nnz; ++k)
      {
        values[k] = dfdy[k] + clock_->alpha * dfdyp[k];
      }
      return 0;
    }

    /**
     * @brief Size the internal variables, the terminals, the inputs and the outputs.
     *
     * @param[in] variables Size of each internal variable, in enumeration order.
     * @param[in] terminals Number of external variables of each terminal, in enumeration order.
     * @param[in] inputs Number of external variables of each input, in enumeration order.
     * @param[in] outputs Number of entries of each output, in enumeration order.
     *
     * @throws std::logic_error if the layout is already set.
     * @throws std::bad_alloc if the local vectors cannot be allocated.
     *
     * @post No terminal, input or output is connected.
     */
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::setLayout(const std::vector<IdxT>& variables,
                                                       const std::vector<IdxT>& terminals,
                                                       const std::vector<IdxT>& inputs,
                                                       const std::vector<IdxT>& outputs)
    {
      // The parameters of a model fix its variables, so a layout is set once
      if (!offsets_.empty())
      {
        throw std::logic_error("EMT::Component: the layout is already set");
      }

      offsets_.assign(1, 0);
      for (const IdxT variable : variables)
      {
        offsets_.push_back(offsets_.back() + variable);
      }
      size_ = offsets_.back();

      // Sizes the connections of one kind, whose offsets start at `first`
      const auto size = [](const std::vector<IdxT>& sizes, std::vector<Connection>& connections, std::vector<IdxT>& offsets, IdxT first)
      {
        IdxT total = 0;
        for (const IdxT entries : sizes)
        {
          connections.push_back({nullptr, std::vector<IdxT>(entries)});
          offsets.push_back(first + total);
          total += entries;
        }
        return total;
      };
      terminals_size_          = size(terminals, terminals_, terminal_offsets_, 0);
      const IdxT inputs_size   = size(inputs, inputs_, input_offsets_, terminals_size_);
      const IdxT outputs_size  = size(outputs, outputs_, output_offsets_, 0);
      const IdxT external_size = terminals_size_ + inputs_size;

      residual_indices_.assign(size_ + terminals_size_ + outputs_size, INVALID_INDEX<IdxT>);
      output_tol_.assign(outputs_size, ZERO<ScalarT>);
      output_values_.assign(outputs_size, ZERO<ScalarT>);
      if (variable_indices_.resize(size_ + external_size) != 0
          || local_y_.resize(external_size) != 0
          || local_yp_.resize(external_size) != 0
          || local_f_.resize(residual_indices_.size()) != 0)
      {
        Log::error() << "EMT::Component: failed to allocate the local variables\n";
        throw std::bad_alloc{};
      }
    }

    /**
     * @brief Local index of external variable `k` of `terminal`.
     *
     * @tparam TerminalT Terminal enumeration of the model.
     * @return Index into the local values and the local residuals.
     */
    template <typename scalar_type, typename index_type>
    template <Utilities::SizedEnum TerminalT>
    typename Component<scalar_type, index_type>::IdxT Component<scalar_type, index_type>::local(TerminalT terminal, IdxT k) const
    {
      return size_ + terminalOffset(terminal) + k;
    }

    /**
     * @brief Construct submodel `symbol`, owned by this component and evaluated by the system.
     *
     * The submodel is named by its symbol, and the system prefixes the id of
     * its parent.
     *
     * @tparam ModelT Model of the submodel.
     * @tparam SymbolT Submodel enumeration of the parent.
     * @tparam Args Constructor argument types of `ModelT`.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT, Utilities::SizedEnum SymbolT, typename... Args>
    ModelT& Component<scalar_type, index_type>::addSubmodel(SymbolT symbol, Args&&... args)
    {
      const auto index = static_cast<IdxT>(symbol);
      if (submodels_.size() <= index)
      {
        submodels_.resize(index + 1);
      }

      auto submodel     = std::make_unique<ModelT>(std::forward<Args>(args)...);
      submodel->id_     = Utilities::enum_name(symbol);
      ModelT& added     = *submodel;
      submodels_[index] = std::move(submodel);
      return added;
    }

    /**
     * @brief Submodel `symbol` as a `ModelT`.
     *
     * @tparam ModelT Model of the submodel.
     * @tparam SymbolT Submodel enumeration of the parent.
     *
     * @pre `addSubmodel` added `symbol`, if at all, as a `ModelT`.
     *
     * @throws std::invalid_argument if the submodel is absent.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT, Utilities::SizedEnum SymbolT>
    const ModelT& Component<scalar_type, index_type>::submodel(SymbolT symbol) const
    {
      const auto index = static_cast<IdxT>(symbol);
      if (index >= submodels_.size() || !submodels_[index])
      {
        throw std::invalid_argument(std::format("{}: submodel '{}' is absent", id_, Utilities::enum_name(symbol)));
      }
      return static_cast<const ModelT&>(*submodels_[index]);
    }

    /**
     * @brief Submodel `symbol` as a `ModelT`, for a caller that sets its values.
     *
     * @tparam ModelT Model of the submodel.
     * @tparam SymbolT Submodel enumeration of the parent.
     *
     * @throws std::invalid_argument if the submodel is absent.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT, Utilities::SizedEnum SymbolT>
    ModelT& Component<scalar_type, index_type>::submodel(SymbolT symbol)
    {
      return const_cast<ModelT&>(std::as_const(*this).template submodel<ModelT>(symbol));
    }

    /**
     * @brief Write the residuals of `model`: of its variables, at its
     *        terminals, and each of its outputs that is connected.
     *
     * A model provides the kernel of each kind it has: `internalResidual`,
     * `externalResidual` if it has terminals, and `output` if it has outputs.
     *
     * @tparam ModelT Model providing the kernels.
     * @tparam T Scalar type of the kernels.
     * @param[in] model This component as its derived type.
     * @param[in] y Internal values.
     * @param[in] yp Internal derivatives.
     * @param[in] w External values.
     * @param[in] wp External derivatives.
     * @param[out] f Local residuals: internal, then those at the terminals and the outputs.
     *
     * @pre The system has resolved the indices of this component.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT, typename T>
    void Component<scalar_type, index_type>::evaluateKernels(
        const ModelT& model, const T* y, const T* yp, const T* w, const T* wp, T* f) const
    {
      model.internalResidual(y, yp, w, wp, f);
      if constexpr (!std::same_as<typename ModelT::TerminalsT, Utilities::EmptyEnum>)
      {
        model.externalResidual(y, yp, w, wp, f + size_);
      }
      if constexpr (!std::same_as<typename ModelT::OutputsT, Utilities::EmptyEnum>)
      {
        for (const auto output : Utilities::enum_values<typename ModelT::OutputsT>())
        {
          if (outputConnected(output))
          {
            model.output(output, y, yp, w, wp, f + size_ + terminals_size_ + outputOffset(output));
          }
        }
      }
    }

    /**
     * @brief Evaluate the kernels of `model` and add the result to the system residual.
     *
     * @tparam ModelT Model providing the kernels.
     * @param[in] model This component as its derived type.
     * @return 0.
     *
     * @pre The system has bound this component and resolved its indices.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    int Component<scalar_type, index_type>::evaluateResidual(const ModelT& model)
    {
      const auto [w, wp] = gather();
      ScalarT* f         = local_f_.getData();

      local_f_.setToZero();
      evaluateKernels(model, y_.getData(), yp_.getData(), w, wp, f);
      scatter(f);
      return 0;
    }

    /**
     * @brief Evaluate every output of `model` at the current values.
     *
     * The outputs are evaluated whether or not they are connected.
     *
     * @tparam ModelT Model providing `output`.
     * @param[in] model This component as its derived type.
     * @return 0.
     *
     * @pre The system has bound this component and resolved its indices.
     * @post `getOutput` gives the entries of each output.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    int Component<scalar_type, index_type>::evaluateOutputs(const ModelT& model)
    {
      const auto [w, wp] = gather();
      std::ranges::fill(output_values_, ZERO<ScalarT>);
      for (const auto output : Utilities::enum_values<typename ModelT::OutputsT>())
      {
        model.output(output, y_.getData(), yp_.getData(), w, wp, output_values_.data() + outputOffset(output));
      }
      return 0;
    }

    /**
     * @brief Copy the external values and derivatives from the system vectors.
     *
     * The internal values and derivatives are read in place, since a bound
     * component views its block of the system vectors.
     *
     * @return The external values and their derivatives, those of the
     *         terminals and then those of the inputs.
     *
     * @pre The system has bound this component and resolved its indices.
     */
    template <typename scalar_type, typename index_type>
    std::pair<const scalar_type*, const scalar_type*> Component<scalar_type, index_type>::gather()
    {
      ScalarT* w  = local_y_.getData();
      ScalarT* wp = local_yp_.getData();
      for (IdxT k = 0; k < external_indices_.size(); ++k)
      {
        w[k]  = system_y_[external_indices_[k]];
        wp[k] = system_yp_[external_indices_[k]];
      }
      return {w, wp};
    }

    /**
     * @brief Add the local residual to the system residual.
     *
     * The residual of an output that is not connected is not added. Devices
     * at one bus add to the same residuals, so a parallel scatter must
     * synchronize.
     *
     * @param[in] f Residual of each local variable.
     */
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::scatter(const ScalarT* f)
    {
      ScalarT* internal = f_.getData();
      for (IdxT j = 0; j < size_; ++j)
      {
        internal[j] += f[j];
      }
      for (const Added& added : added_)
      {
        system_f_[added.residual] += f[size_ + added.external];
      }
    }

    /**
     * @brief Resolve the system index of each external value and residual.
     *
     * A connection to an external variable of another component resolves
     * through that component's indices.
     *
     * @pre The system has bound every component, and `setIndices` has run on
     *      each component whose external variable a terminal, an input or an
     *      output is connected to.
     * @post `variable_indices_` holds the system variable of every local
     *       value, and `residual_indices_` the system residual of every local
     *       residual, with no index for an output that is not connected.
     */
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::setIndices()
    {
      IdxT* variables = variable_indices_.getData();

      external_indices_.clear();
      added_.clear();
      IdxT value    = size_;
      IdxT residual = size_;
      for (const Connection& terminal : terminals_)
      {
        for (const IdxT i : terminal.indices)
        {
          variables[value]            = terminal.component->variable_indices_.getData()[i];
          residual_indices_[residual] = terminal.component->residual_indices_[i];
          external_indices_.push_back(variables[value]);
          added_.push_back({residual - size_, residual_indices_[residual]});
          ++value;
          ++residual;
        }
      }
      for (const Connection& input : inputs_)
      {
        for (const IdxT i : input.indices)
        {
          variables[value] = input.component->variable_indices_.getData()[i];
          external_indices_.push_back(variables[value]);
          ++value;
        }
      }
      for (const Connection& output : outputs_)
      {
        for (const IdxT i : output.indices)
        {
          if (output.component != nullptr)
          {
            residual_indices_[residual] = output.component->residual_indices_[i];
            added_.push_back({residual - size_, residual_indices_[residual]});
          }
          ++residual;
        }
      }
    }
  } // namespace EMT
} // namespace GridKit

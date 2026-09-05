#pragma once

#include <numeric>

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobian.hpp>
#include <GridKit/Model/EMT/Component.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Compose an external residual derivative with its signal gradients.
     *
     * Enzyme differentiates the local residual kernel. Computed signals expand
     * each input column into its global DAE columns outside that kernel.
     * Callers reserve their usual dense capacity times externalJacobianExpansion().
     */
    template <typename ModelT, Enzyme::Sparse::Equation equation, Enzyme::Sparse::Variable variable>
    struct SignalJacobian
    {
      using ScalarT    = typename ModelT::ScalarT;
      using IdxT       = typename ModelT::IdxT;
      using RealT      = typename ModelT::RealT;
      using ComponentT = Component<ScalarT, IdxT>;
      using SignalT    = Signal<ScalarT, IdxT>;

      static void eval(ComponentT* owner, ModelT* model, size_t n_res, size_t n_var, const IdxT* res_indices, const IdxT* var_indices, const ScalarT* y, const ScalarT* yp, const ScalarT* y_ext, const ScalarT* yp_ext, IdxT* rows, IdxT* cols, RealT* vals, IdxT& nnz, RealT scaling = RealT{1})
      {
        using Kernel = Enzyme::Sparse::SparseJacobian<ModelT, equation, variable>;
        if (!owner->hasComputedInputs())
        {
          Kernel::eval(model, n_res, n_var, res_indices, var_indices, y, yp, y_ext, yp_ext, rows, cols, vals, nnz, scaling);
          return;
        }
        std::vector<IdxT> local_indices(n_var);
        std::iota(local_indices.begin(), local_indices.end(), IdxT{0});
        std::vector<IdxT>  local_rows(n_res * n_var), local_cols(n_res * n_var);
        std::vector<RealT> local_vals(n_res * n_var);
        IdxT               local_nnz = 0;
        Kernel::eval(model, n_res, n_var, res_indices, local_indices.data(), y, yp, y_ext, yp_ext, local_rows.data(), local_cols.data(), local_vals.data(), local_nnz, scaling);

        std::vector<typename SignalT::GradientT> gradients(n_var);
        const auto&                              signals = owner->externalVariableSignals();
        for (size_t n = 0; n < n_var; ++n)
        {
          if (signals[n] != nullptr && signals[n]->computed())
          {
            if constexpr (variable == Enzyme::Sparse::Variable::YpExt)
            {
              // Computed algebraic signals expose no derivative input.
              continue;
            }
            signals[n]->appendGradient(gradients[n]);
          }
          else if (var_indices[n] != INVALID_INDEX<IdxT>)
          {
            gradients[n].emplace_back(var_indices[n], RealT{1});
          }
        }
        for (size_t j = 0; j < static_cast<size_t>(local_nnz); ++j)
        {
          for (const auto& [column, coefficient] : gradients[static_cast<size_t>(local_cols[j])])
          {
            rows[nnz] = local_rows[j];
            cols[nnz] = column;
            vals[nnz] = local_vals[j] * coefficient;
            ++nnz;
          }
        }
      }
    };
  } // namespace EMT
} // namespace GridKit

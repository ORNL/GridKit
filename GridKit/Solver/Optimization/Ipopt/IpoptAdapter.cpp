#include "IpoptAdapter.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace
{
  template <typename target_type, typename value_type>
  constexpr bool isNonnegativeAndInRange(value_type value)
  {
    if constexpr (std::is_signed_v<value_type>)
    {
      if (value < 0)
      {
        return false;
      }
    }
    return std::in_range<target_type>(value);
  }
} // namespace

namespace GridKit
{
  namespace Optimization
  {
    template <class scalar_type, typename index_type>
    IpoptAdapter<scalar_type, index_type>::IpoptAdapter(EvaluatorT& model)
      : model_(model)
    {
      static_assert(std::is_same_v<ScalarT, Ipopt::Number>,
                    "The Ipopt adapter currently supports double variables only");
      static_assert(std::is_same_v<RealT, Ipopt::Number>,
                    "The Ipopt adapter currently supports double values only");
      static_assert(std::is_integral_v<IdxT>,
                    "The Ipopt adapter requires an integral index type");

      if (!model_.hasJacobian())
      {
        throw std::invalid_argument("IpoptAdapter requires an exact sparse Jacobian");
      }
      if (!model_.hasHessian())
      {
        throw std::invalid_argument("IpoptAdapter requires an exact sparse Hessian");
      }

      auto* jacobian = model_.getCsrJacobian();
      auto* hessian  = model_.getCsrHessian();
      if (jacobian == nullptr || hessian == nullptr)
      {
        throw std::invalid_argument("IpoptAdapter requires allocated derivative structures");
      }

      const IdxT model_variable_count    = model_.size();
      const IdxT model_constraint_count  = model_.sizeConstraints();
      const IdxT model_jacobian_nonzeros = jacobian->getNnz();
      const IdxT model_hessian_nonzeros  = hessian->getNnz();

      if (!isNonnegativeAndInRange<Index>(model_variable_count)
          || !isNonnegativeAndInRange<Index>(model_constraint_count)
          || !isNonnegativeAndInRange<Index>(model_jacobian_nonzeros)
          || !isNonnegativeAndInRange<Index>(model_hessian_nonzeros))
      {
        throw std::invalid_argument("IpoptAdapter model dimensions exceed Ipopt index limits");
      }

      variable_count_    = static_cast<Index>(model_variable_count);
      constraint_count_  = static_cast<Index>(model_constraint_count);
      jacobian_nonzeros_ = static_cast<Index>(model_jacobian_nonzeros);
      hessian_nonzeros_  = static_cast<Index>(model_hessian_nonzeros);

      if (model_.variables().getSize() != model_variable_count
          || model_.variableLowerBounds().getSize() != model_variable_count
          || model_.variableUpperBounds().getSize() != model_variable_count
          || model_.objectiveGradient().getSize() != model_variable_count
          || model_.constraints().getSize() != model_constraint_count
          || model_.constraintLowerBounds().getSize() != model_constraint_count
          || model_.constraintUpperBounds().getSize() != model_constraint_count)
      {
        throw std::invalid_argument("IpoptAdapter model vector dimensions are inconsistent");
      }

      if (jacobian->getNumRows() != model_constraint_count
          || jacobian->getNumColumns() != model_variable_count
          || hessian->getNumRows() != model_variable_count
          || hessian->getNumColumns() != model_variable_count)
      {
        throw std::invalid_argument("IpoptAdapter derivative dimensions do not match the model");
      }

      const auto validate_structure =
          [](CsrMatrixT* matrix, bool lower_triangular)
      {
        const IdxT  row_count     = matrix->getNumRows();
        const IdxT  column_count  = matrix->getNumColumns();
        const IdxT  nonzero_count = matrix->getNnz();
        const auto* rows          = matrix->getRowData(memory::HOST);
        const auto* columns       = matrix->getColData(memory::HOST);
        const auto* values        = matrix->getValues(memory::HOST);

        if (rows == nullptr
            || (nonzero_count > 0
                && (columns == nullptr || values == nullptr))
            || rows[0] != 0)
        {
          return false;
        }

        for (std::size_t row = 0;
             row < static_cast<std::size_t>(row_count);
             ++row)
        {
          const IdxT begin = rows[row];
          const IdxT end   = rows[row + 1];
          if (!isNonnegativeAndInRange<std::size_t>(begin)
              || !isNonnegativeAndInRange<std::size_t>(end)
              || begin > end
              || end > nonzero_count)
          {
            return false;
          }

          for (std::size_t entry = static_cast<std::size_t>(begin);
               entry < static_cast<std::size_t>(end);
               ++entry)
          {
            const IdxT column = columns[entry];
            if (!isNonnegativeAndInRange<Index>(column)
                || column >= column_count
                || (entry > static_cast<std::size_t>(begin)
                    && column <= columns[entry - 1])
                || (lower_triangular
                    && column > static_cast<IdxT>(row)))
            {
              return false;
            }
          }
        }
        return rows[static_cast<std::size_t>(row_count)] == nonzero_count;
      };

      if (!validate_structure(jacobian, false))
      {
        throw std::invalid_argument("IpoptAdapter requires a valid deduplicated CSR Jacobian");
      }
      if (!validate_structure(hessian, true))
      {
        throw std::invalid_argument(
            "IpoptAdapter requires a valid lower-triangular CSR Hessian");
      }
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::get_nlp_info(
        Index&          variable_count,
        Index&          constraint_count,
        Index&          jacobian_nonzeros,
        Index&          hessian_nonzeros,
        IndexStyleEnum& index_style)
    {
      variable_count    = variable_count_;
      constraint_count  = constraint_count_;
      jacobian_nonzeros = jacobian_nonzeros_;
      hessian_nonzeros  = hessian_nonzeros_;
      index_style       = C_STYLE;
      return true;
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::get_bounds_info(
        Index   variable_count,
        Number* variable_lower_bounds,
        Number* variable_upper_bounds,
        Index   constraint_count,
        Number* constraint_lower_bounds,
        Number* constraint_upper_bounds)
    {
      if (variable_count != variable_count_
          || constraint_count != constraint_count_)
      {
        return false;
      }

      std::copy_n(model_.variableLowerBounds().getData(memory::HOST),
                  static_cast<std::size_t>(variable_count),
                  variable_lower_bounds);
      std::copy_n(model_.variableUpperBounds().getData(memory::HOST),
                  static_cast<std::size_t>(variable_count),
                  variable_upper_bounds);
      std::copy_n(model_.constraintLowerBounds().getData(memory::HOST),
                  static_cast<std::size_t>(constraint_count),
                  constraint_lower_bounds);
      std::copy_n(model_.constraintUpperBounds().getData(memory::HOST),
                  static_cast<std::size_t>(constraint_count),
                  constraint_upper_bounds);
      return true;
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::get_starting_point(
        Index                    variable_count,
        bool                     initialize_variables,
        Number*                  variables,
        bool                     initialize_bound_multipliers,
        [[maybe_unused]] Number* lower_bound_multipliers,
        [[maybe_unused]] Number* upper_bound_multipliers,
        Index                    constraint_count,
        bool                     initialize_constraint_multipliers,
        [[maybe_unused]] Number* constraint_multipliers)
    {
      if (variable_count != variable_count_
          || constraint_count != constraint_count_
          || !initialize_variables
          || initialize_bound_multipliers
          || initialize_constraint_multipliers)
      {
        return false;
      }

      std::copy_n(model_.variables().getData(memory::HOST),
                  static_cast<std::size_t>(variable_count),
                  variables);
      return true;
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::eval_f(
        Index                 variable_count,
        const Number*         variables,
        [[maybe_unused]] bool new_variables,
        Number&               objective)
    {
      if (!syncVariables(variable_count, variables)
          || model_.evaluateObjective() != 0)
      {
        return false;
      }
      objective = model_.objective();
      return true;
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::eval_grad_f(
        Index                 variable_count,
        const Number*         variables,
        [[maybe_unused]] bool new_variables,
        Number*               gradient)
    {
      if (!syncVariables(variable_count, variables)
          || model_.evaluateObjectiveGradient() != 0)
      {
        return false;
      }

      std::copy_n(model_.objectiveGradient().getData(memory::HOST),
                  static_cast<std::size_t>(variable_count),
                  gradient);
      return true;
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::eval_g(
        Index                 variable_count,
        const Number*         variables,
        [[maybe_unused]] bool new_variables,
        Index                 constraint_count,
        Number*               constraints)
    {
      if (constraint_count != constraint_count_
          || !syncVariables(variable_count, variables)
          || model_.evaluateConstraints() != 0)
      {
        return false;
      }

      std::copy_n(model_.constraints().getData(memory::HOST),
                  static_cast<std::size_t>(constraint_count),
                  constraints);
      return true;
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::eval_jac_g(
        Index                 variable_count,
        const Number*         variables,
        [[maybe_unused]] bool new_variables,
        Index                 constraint_count,
        Index                 nonzero_count,
        Index*                rows,
        Index*                columns,
        Number*               values)
    {
      auto* jacobian = model_.getCsrJacobian();
      if (jacobian == nullptr
          || variable_count != variable_count_
          || constraint_count != constraint_count_
          || nonzero_count != jacobian_nonzeros_
          || jacobian->getNnz() != static_cast<IdxT>(jacobian_nonzeros_))
      {
        return false;
      }

      if (values == nullptr)
      {
        if (nonzero_count > 0 && (rows == nullptr || columns == nullptr))
        {
          return false;
        }

        const auto* row_data    = jacobian->getRowData(memory::HOST);
        const auto* column_data = jacobian->getColData(memory::HOST);
        Index       output      = 0;
        for (IdxT row = 0; row < jacobian->getNumRows(); ++row)
        {
          for (IdxT entry = row_data[static_cast<std::size_t>(row)];
               entry < row_data[static_cast<std::size_t>(row) + 1];
               ++entry)
          {
            rows[static_cast<std::size_t>(output)] = static_cast<Index>(row);
            columns[static_cast<std::size_t>(output)] =
                static_cast<Index>(column_data[static_cast<std::size_t>(entry)]);
            ++output;
          }
        }
        return output == nonzero_count;
      }

      if (!syncVariables(variable_count, variables)
          || model_.evaluateJacobian() != 0)
      {
        return false;
      }

      std::copy_n(jacobian->getValues(memory::HOST),
                  static_cast<std::size_t>(nonzero_count),
                  values);
      return true;
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::eval_h(
        Index                 variable_count,
        const Number*         variables,
        [[maybe_unused]] bool new_variables,
        Number                objective_factor,
        Index                 constraint_count,
        const Number*         constraint_multipliers,
        [[maybe_unused]] bool new_constraint_multipliers,
        Index                 nonzero_count,
        Index*                rows,
        Index*                columns,
        Number*               values)
    {
      auto* hessian = model_.getCsrHessian();
      if (hessian == nullptr
          || variable_count != variable_count_
          || constraint_count != constraint_count_
          || nonzero_count != hessian_nonzeros_
          || hessian->getNnz() != static_cast<IdxT>(hessian_nonzeros_))
      {
        return false;
      }

      if (values == nullptr)
      {
        if (nonzero_count > 0 && (rows == nullptr || columns == nullptr))
        {
          return false;
        }

        const auto* row_data    = hessian->getRowData(memory::HOST);
        const auto* column_data = hessian->getColData(memory::HOST);
        Index       output      = 0;
        for (IdxT row = 0; row < hessian->getNumRows(); ++row)
        {
          for (IdxT entry = row_data[static_cast<std::size_t>(row)];
               entry < row_data[static_cast<std::size_t>(row) + 1];
               ++entry)
          {
            rows[static_cast<std::size_t>(output)] = static_cast<Index>(row);
            columns[static_cast<std::size_t>(output)] =
                static_cast<Index>(column_data[static_cast<std::size_t>(entry)]);
            ++output;
          }
        }
        return output == nonzero_count;
      }

      if (constraint_count > 0 && constraint_multipliers == nullptr)
      {
        return false;
      }
      if (!syncVariables(variable_count, variables)
          || model_.evaluateHessian(
                 objective_factor,
                 constraint_multipliers,
                 static_cast<IdxT>(constraint_count))
                 != 0)
      {
        return false;
      }

      std::copy_n(hessian->getValues(memory::HOST),
                  static_cast<std::size_t>(nonzero_count),
                  values);
      return true;
    }

    template <class scalar_type, typename index_type>
    void IpoptAdapter<scalar_type, index_type>::finalize_solution(
        [[maybe_unused]] SolverReturn               status,
        Index                                       variable_count,
        const Number*                               variables,
        [[maybe_unused]] const Number*              lower_bound_multipliers,
        [[maybe_unused]] const Number*              upper_bound_multipliers,
        [[maybe_unused]] Index                      constraint_count,
        [[maybe_unused]] const Number*              constraints,
        [[maybe_unused]] const Number*              constraint_multipliers,
        [[maybe_unused]] Number                     objective,
        [[maybe_unused]] const IpoptData*           ipopt_data,
        [[maybe_unused]] IpoptCalculatedQuantities* ipopt_quantities)
    {
      syncVariables(variable_count, variables);
    }

    template <class scalar_type, typename index_type>
    bool IpoptAdapter<scalar_type, index_type>::syncVariables(
        Index variable_count, const Number* variables)
    {
      if (variables == nullptr
          || variable_count != variable_count_)
      {
        return false;
      }
      return model_.variables().copyFromExternal(variables, memory::HOST, memory::HOST) == 0;
    }

    template class IpoptAdapter<Ipopt::Number, int>;
    template class IpoptAdapter<Ipopt::Number, long int>;
    template class IpoptAdapter<Ipopt::Number, std::size_t>;
  } // namespace Optimization
} // namespace GridKit

#include "IpoptSolver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

#include <IpIpoptApplication.hpp>
#include <IpIpoptData.hpp>
#include <IpTNLP.hpp>

namespace AnalysisManager
{
  namespace IpoptInterface
  {
    namespace
    {
      constexpr Ipopt::Number IPOPT_LOWER_INFINITY = -1.0e19;
      constexpr Ipopt::Number IPOPT_UPPER_INFINITY = 1.0e19;

      OptimizationStatus translateStatus(Ipopt::SolverReturn status)
      {
        switch (status)
        {
        case Ipopt::SUCCESS:
          return OptimizationStatus::SOLVED;
        case Ipopt::STOP_AT_ACCEPTABLE_POINT:
        case Ipopt::FEASIBLE_POINT_FOUND:
          return OptimizationStatus::ACCEPTABLE;
        case Ipopt::LOCAL_INFEASIBILITY:
          return OptimizationStatus::INFEASIBLE;
        case Ipopt::MAXITER_EXCEEDED:
        case Ipopt::CPUTIME_EXCEEDED:
        case Ipopt::WALLTIME_EXCEEDED:
          return OptimizationStatus::LIMIT_REACHED;
        case Ipopt::USER_REQUESTED_STOP:
          return OptimizationStatus::USER_STOPPED;
        case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM:
        case Ipopt::INVALID_OPTION:
          return OptimizationStatus::INVALID_PROBLEM;
        default:
          return OptimizationStatus::SOLVER_ERROR;
        }
      }

      OptimizationStatus translateStatus(Ipopt::ApplicationReturnStatus status)
      {
        switch (status)
        {
        case Ipopt::Solve_Succeeded:
          return OptimizationStatus::SOLVED;
        case Ipopt::Solved_To_Acceptable_Level:
        case Ipopt::Feasible_Point_Found:
          return OptimizationStatus::ACCEPTABLE;
        case Ipopt::Infeasible_Problem_Detected:
          return OptimizationStatus::INFEASIBLE;
        case Ipopt::Maximum_Iterations_Exceeded:
        case Ipopt::Maximum_CpuTime_Exceeded:
        case Ipopt::Maximum_WallTime_Exceeded:
          return OptimizationStatus::LIMIT_REACHED;
        case Ipopt::User_Requested_Stop:
          return OptimizationStatus::USER_STOPPED;
        case Ipopt::Not_Enough_Degrees_Of_Freedom:
        case Ipopt::Invalid_Problem_Definition:
        case Ipopt::Invalid_Option:
          return OptimizationStatus::INVALID_PROBLEM;
        default:
          return OptimizationStatus::SOLVER_ERROR;
        }
      }
    } // namespace

    template <class ScalarT, typename IdxT>
    class Solver<ScalarT, IdxT>::Callbacks : public Ipopt::TNLP
    {
      using RealT       = typename Solver::EvaluatorT::RealT;
      using RealVectorT = typename Solver::EvaluatorT::RealVectorT;
      using VectorT     = typename Solver::EvaluatorT::VectorT;
      using CsrMatrixT  = typename Solver::EvaluatorT::CsrMatrixT;

      static_assert(std::is_same_v<RealT, Ipopt::Number>,
                    "Ipopt requires the evaluator real type to match Ipopt::Number");

    public:
      Callbacks(EvaluatorT* model, ResultT& result)
        : model_(model),
          result_(result)
      {
      }

      bool get_nlp_info(Ipopt::Index&   n,
                        Ipopt::Index&   m,
                        Ipopt::Index&   nnz_jac_g,
                        Ipopt::Index&   nnz_h_lag,
                        IndexStyleEnum& index_style) override
      {
        if (!toIpoptIndex(model_->size(), n) || !toIpoptIndex(model_->sizeResidual(), m) || !toIpoptIndex(model_->nnz(), nnz_jac_g))
        {
          return false;
        }

        nnz_h_lag   = 0;
        index_style = C_STYLE;
        if (n <= 0 || (m == 0 && nnz_jac_g != 0))
        {
          return false;
        }
        return m == 0 || model_->hasJacobian();
      }

      bool get_bounds_info(Ipopt::Index   n,
                           Ipopt::Number* x_l,
                           Ipopt::Number* x_u,
                           Ipopt::Index   m,
                           Ipopt::Number* g_l,
                           Ipopt::Number* g_u) override
      {
        if (!matchesDimensions(n, m) || model_->getVariableBounds(variable_lower_, variable_upper_) != 0 || model_->getResidualBounds(residual_lower_, residual_upper_) != 0 || variable_lower_.getSize() != model_->size() || variable_upper_.getSize() != model_->size() || residual_lower_.getSize() != model_->sizeResidual() || residual_upper_.getSize() != model_->sizeResidual())
        {
          return false;
        }

        const auto* variable_lower = variable_lower_.getData();
        const auto* variable_upper = variable_upper_.getData();
        for (Ipopt::Index i = 0; i < n; ++i)
        {
          x_l[i] = lowerBound(variable_lower[i]);
          x_u[i] = upperBound(variable_upper[i]);
        }

        const auto* residual_lower = residual_lower_.getData();
        const auto* residual_upper = residual_upper_.getData();
        for (Ipopt::Index i = 0; i < m; ++i)
        {
          g_l[i] = lowerBound(residual_lower[i]);
          g_u[i] = upperBound(residual_upper[i]);
        }
        return true;
      }

      bool get_starting_point(Ipopt::Index                    n,
                              bool                            init_x,
                              Ipopt::Number*                  x,
                              bool                            init_z,
                              [[maybe_unused]] Ipopt::Number* z_L,
                              [[maybe_unused]] Ipopt::Number* z_U,
                              Ipopt::Index                    m,
                              bool                            init_lambda,
                              [[maybe_unused]] Ipopt::Number* lambda) override
      {
        if (!matchesDimensions(n, m) || !init_x || init_z || init_lambda || model_->y().getSize() != model_->size())
        {
          return false;
        }

        const auto* variables = model_->y().getData();
        for (Ipopt::Index i = 0; i < n; ++i)
        {
          x[i] = static_cast<Ipopt::Number>(variables[i]);
        }
        return true;
      }

      bool eval_f(Ipopt::Index         n,
                  const Ipopt::Number* x,
                  bool                 new_x,
                  Ipopt::Number&       obj_value) override
      {
        if (!matchesVariables(n) || !syncVariables(x, new_x))
        {
          return false;
        }

        if (!model_->hasObjective())
        {
          obj_value = 0.0;
          return true;
        }

        if (!objective_current_)
        {
          if (model_->evaluateObjective() != 0)
          {
            return false;
          }
          objective_current_ = true;
        }

        obj_value = static_cast<Ipopt::Number>(model_->objective());
        return true;
      }

      bool eval_grad_f(Ipopt::Index         n,
                       const Ipopt::Number* x,
                       bool                 new_x,
                       Ipopt::Number*       grad_f) override
      {
        if (!matchesVariables(n) || !syncVariables(x, new_x))
        {
          return false;
        }

        if (!model_->hasObjective())
        {
          std::fill_n(grad_f, n, 0.0);
          return true;
        }

        if (!evaluateDerivatives())
        {
          return false;
        }

        const VectorT& gradient = model_->getObjectiveGradient();
        if (gradient.getSize() != model_->size())
        {
          return false;
        }

        const auto* values = gradient.getData();
        for (Ipopt::Index i = 0; i < n; ++i)
        {
          grad_f[i] = static_cast<Ipopt::Number>(values[i]);
        }
        return true;
      }

      bool eval_g(Ipopt::Index         n,
                  const Ipopt::Number* x,
                  bool                 new_x,
                  Ipopt::Index         m,
                  Ipopt::Number*       g) override
      {
        if (!matchesDimensions(n, m))
        {
          return false;
        }

        if (!syncVariables(x, new_x))
        {
          return false;
        }

        if (m == 0)
        {
          return true;
        }

        if (!residual_current_)
        {
          if (model_->evaluateResidual() != 0)
          {
            return false;
          }
          residual_current_ = true;
        }

        const VectorT& residual = model_->getResidual();
        if (residual.getSize() != model_->sizeResidual())
        {
          return false;
        }

        const auto* values = residual.getData();
        for (Ipopt::Index i = 0; i < m; ++i)
        {
          g[i] = static_cast<Ipopt::Number>(values[i]);
        }
        return true;
      }

      bool eval_jac_g(Ipopt::Index         n,
                      const Ipopt::Number* x,
                      bool                 new_x,
                      Ipopt::Index         m,
                      Ipopt::Index         nele_jac,
                      Ipopt::Index*        iRow,
                      Ipopt::Index*        jCol,
                      Ipopt::Number*       values) override
      {
        if (!matchesDimensions(n, m) || !matchesNonzeros(nele_jac))
        {
          return false;
        }

        if (m == 0 && nele_jac == 0)
        {
          return x == nullptr || syncVariables(x, new_x);
        }

        CsrMatrixT* jacobian = model_->getCsrJacobian();
        if (!validJacobian(jacobian))
        {
          return false;
        }

        if (values == nullptr)
        {
          const auto* rows = jacobian->getRowData();
          const auto* cols = jacobian->getColData();
          for (Ipopt::Index row = 0; row < m; ++row)
          {
            const IdxT begin = rows[row];
            const IdxT end   = rows[row + 1];
            for (IdxT entry = begin; entry < end; ++entry)
            {
              Ipopt::Index column{};
              Ipopt::Index slot{};
              if (!toIpoptIndex(cols[entry], column) || !toIpoptIndex(entry, slot))
              {
                return false;
              }
              iRow[slot] = row;
              jCol[slot] = column;
            }
          }
          return true;
        }

        if (!syncVariables(x, new_x) || !evaluateDerivatives())
        {
          return false;
        }

        const auto* jacobian_values = jacobian->getValues();
        for (Ipopt::Index i = 0; i < nele_jac; ++i)
        {
          values[i] = static_cast<Ipopt::Number>(jacobian_values[i]);
        }
        return true;
      }

      bool eval_h([[maybe_unused]] Ipopt::Index         n,
                  [[maybe_unused]] const Ipopt::Number* x,
                  [[maybe_unused]] bool                 new_x,
                  [[maybe_unused]] Ipopt::Number        obj_factor,
                  [[maybe_unused]] Ipopt::Index         m,
                  [[maybe_unused]] const Ipopt::Number* lambda,
                  [[maybe_unused]] bool                 new_lambda,
                  [[maybe_unused]] Ipopt::Index         nele_hess,
                  [[maybe_unused]] Ipopt::Index*        iRow,
                  [[maybe_unused]] Ipopt::Index*        jCol,
                  [[maybe_unused]] Ipopt::Number*       values) override
      {
        return false;
      }

      void finalize_solution(Ipopt::SolverReturn                                status,
                             Ipopt::Index                                       n,
                             const Ipopt::Number*                               x,
                             const Ipopt::Number*                               z_L,
                             const Ipopt::Number*                               z_U,
                             Ipopt::Index                                       m,
                             const Ipopt::Number*                               g,
                             const Ipopt::Number*                               lambda,
                             Ipopt::Number                                      obj_value,
                             const Ipopt::IpoptData*                            ip_data,
                             [[maybe_unused]] Ipopt::IpoptCalculatedQuantities* ip_cq) override
      {
        result_.status    = translateStatus(status);
        result_.objective = static_cast<RealT>(obj_value);

        if (ip_data != nullptr)
        {
          result_.iterations = static_cast<IdxT>(ip_data->iter_count());
        }

        if (n > 0)
        {
          result_.primal.assign(x, x + n);
          result_.lower_bound_duals.assign(z_L, z_L + n);
          result_.upper_bound_duals.assign(z_U, z_U + n);
          if (!syncVariables(x, true))
          {
            result_.status = OptimizationStatus::SOLVER_ERROR;
          }
        }

        if (m > 0)
        {
          result_.constraints.assign(g, g + m);
          result_.constraint_duals.assign(lambda, lambda + m);
          setConstraintViolation(g, m);
        }
      }

    private:
      template <typename SourceIndexT>
      bool toIpoptIndex(SourceIndexT source, Ipopt::Index& destination) const
      {
        destination = static_cast<Ipopt::Index>(source);
        if (destination < 0)
        {
          return false;
        }
        return static_cast<SourceIndexT>(destination) == source;
      }

      bool matchesVariables(Ipopt::Index n) const
      {
        Ipopt::Index expected{};
        return toIpoptIndex(model_->size(), expected) && n == expected;
      }

      bool matchesDimensions(Ipopt::Index n, Ipopt::Index m) const
      {
        Ipopt::Index expected_residuals{};
        return matchesVariables(n) && toIpoptIndex(model_->sizeResidual(), expected_residuals) && m == expected_residuals;
      }

      bool matchesNonzeros(Ipopt::Index nnz) const
      {
        Ipopt::Index expected{};
        return toIpoptIndex(model_->nnz(), expected) && nnz == expected;
      }

      bool syncVariables(const Ipopt::Number* x, bool new_x)
      {
        if (!new_x && variables_current_)
        {
          return true;
        }

        VectorT& variables = model_->y();
        if (variables.getSize() != model_->size())
        {
          return false;
        }

        auto* values = variables.getData();
        for (IdxT i = 0; i < model_->size(); ++i)
        {
          values[i] = static_cast<ScalarT>(x[i]);
        }
        if (variables.setDataUpdated() != 0)
        {
          return false;
        }

        variables_current_   = true;
        objective_current_   = false;
        residual_current_    = false;
        derivatives_current_ = false;
        return true;
      }

      bool evaluateDerivatives()
      {
        if (derivatives_current_)
        {
          return true;
        }

        if (model_->evaluateJacobian() != 0)
        {
          return false;
        }
        derivatives_current_ = true;
        return true;
      }

      bool validJacobian(CsrMatrixT* jacobian) const
      {
        if (jacobian == nullptr || jacobian->getNumRows() != model_->sizeResidual() || jacobian->getNumColumns() != model_->size() || jacobian->getNnz() != model_->nnz())
        {
          return false;
        }

        const auto* rows   = jacobian->getRowData();
        const auto* cols   = jacobian->getColData();
        const auto* values = jacobian->getValues();
        if (rows == nullptr || (model_->nnz() > 0 && (cols == nullptr || values == nullptr)))
        {
          return false;
        }

        if (rows[0] != 0 || rows[model_->sizeResidual()] != model_->nnz())
        {
          return false;
        }

        for (IdxT row = 0; row < model_->sizeResidual(); ++row)
        {
          const IdxT begin = rows[row];
          const IdxT end   = rows[row + 1];
          if (isNegative(begin) || isNegative(end) || begin > end || end > model_->nnz())
          {
            return false;
          }
        }

        for (IdxT entry = 0; entry < model_->nnz(); ++entry)
        {
          if (isNegative(cols[entry]) || cols[entry] >= model_->size())
          {
            return false;
          }
        }
        return true;
      }

      static bool isNegative(IdxT value)
      {
        if constexpr (std::is_signed_v<IdxT>)
        {
          return value < 0;
        }
        return false;
      }

      static Ipopt::Number lowerBound(RealT value)
      {
        return std::isinf(value) && value < 0.0
                   ? IPOPT_LOWER_INFINITY
                   : static_cast<Ipopt::Number>(value);
      }

      static Ipopt::Number upperBound(RealT value)
      {
        return std::isinf(value) && value > 0.0
                   ? IPOPT_UPPER_INFINITY
                   : static_cast<Ipopt::Number>(value);
      }

      void setConstraintViolation(const Ipopt::Number* values, Ipopt::Index m)
      {
        if (residual_lower_.getSize() != model_->sizeResidual() || residual_upper_.getSize() != model_->sizeResidual())
        {
          return;
        }

        const auto* lower = residual_lower_.getData();
        const auto* upper = residual_upper_.getData();
        RealT       violation{};
        for (Ipopt::Index i = 0; i < m; ++i)
        {
          violation = std::max(violation, lower[i] - values[i]);
          violation = std::max(violation, values[i] - upper[i]);
        }
        result_.constraint_violation = violation;
      }

    private:
      EvaluatorT* model_{nullptr};
      ResultT&    result_;

      RealVectorT variable_lower_;
      RealVectorT variable_upper_;
      RealVectorT residual_lower_;
      RealVectorT residual_upper_;

      bool variables_current_{false};
      bool objective_current_{false};
      bool residual_current_{false};
      bool derivatives_current_{false};
    };

    template <class ScalarT, typename IdxT>
    Solver<ScalarT, IdxT>::Solver(EvaluatorT* model)
      : OptimizationSolver<ScalarT, IdxT>(model)
    {
      if (model == nullptr)
      {
        throw std::invalid_argument("Ipopt solver requires a model evaluator");
      }
    }

    template <class ScalarT, typename IdxT>
    typename Solver<ScalarT, IdxT>::ResultT Solver<ScalarT, IdxT>::solve(const OptionsT& options)
    {
      ResultT result;

      Ipopt::SmartPtr<Callbacks> callbacks   = new Callbacks(model_, result);
      auto                       application = IpoptApplicationFactory();

      bool options_valid  = application->Options()->SetNumericValue("tol", options.tolerance);
      options_valid      &= application->Options()->SetIntegerValue("max_iter", options.maximum_iterations);
      options_valid      &= application->Options()->SetIntegerValue("print_level", options.print_level);
      options_valid      &= application->Options()->SetStringValue("hessian_approximation", "limited-memory");
      if (!options.linear_solver.empty())
      {
        options_valid &= application->Options()->SetStringValue("linear_solver", options.linear_solver);
      }
      if (!options_valid)
      {
        result.status = OptimizationStatus::INVALID_PROBLEM;
        return result;
      }

      const auto initialize_status = application->Initialize();
      if (initialize_status != Ipopt::Solve_Succeeded)
      {
        result.status = translateStatus(initialize_status);
        return result;
      }

      const auto solve_status = application->OptimizeTNLP(callbacks);
      if (result.status == OptimizationStatus::NOT_RUN)
      {
        result.status = translateStatus(solve_status);
      }
      return result;
    }

    template class Solver<Ipopt::Number, int>;
    template class Solver<Ipopt::Number, long int>;
    template class Solver<Ipopt::Number, std::size_t>;

  } // namespace IpoptInterface
} // namespace AnalysisManager

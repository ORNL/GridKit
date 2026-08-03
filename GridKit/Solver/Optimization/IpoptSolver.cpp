/**
 * @file IpoptSolver.cpp
 *
 * @brief Interior-point solver backend wrapping Ipopt.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "IpoptSolver.hpp"

#include <sstream>
#include <type_traits>

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpTNLP.hpp>

namespace
{
  /**
   * @brief Presents an OptimizationModel to Ipopt as a TNLP.
   *
   * The dense row-major Jacobian maps one-to-one onto Ipopt's triplet
   * structure enumerated row-major, and the dense symmetric Hessian is
   * exposed through its lower triangle. The solution vector doubles as
   * the starting point; finalize_solution writes the result back into it.
   *
   * Ipopt treats a false evaluator return as a recoverable point
   * failure, so the model's negative (fatal) returns are recorded here:
   * the intermediate callback then stops the iteration and the solver
   * reports the failure through its own error code.
   */
  template <typename scalar_type, typename index_type>
  class ModelAdapter final : public Ipopt::TNLP
  {
  public:
    using ModelT =
        GridKit::Optimization::OptimizationModel<scalar_type, index_type>;
    using StatsT = typename GridKit::Optimization::
        IpoptSolver<scalar_type, index_type>::Stats;
    using RealT = typename ModelT::RealT;

    static_assert(std::is_same<RealT, Ipopt::Number>::value,
                  "The Ipopt backend requires RealT to be Ipopt::Number");

    ModelAdapter(const ModelT&       model,
                 std::vector<RealT>& solution,
                 StatsT&             stats,
                 bool                exact_hessian)
      : model_(model),
        solution_(solution),
        stats_(stats),
        exact_hessian_(exact_hessian),
        hessian_buffer_(static_cast<size_t>(model.variableCount()) * static_cast<size_t>(model.variableCount()),
                        RealT{0})
    {
    }

    bool get_nlp_info(Ipopt::Index&   n,
                      Ipopt::Index&   m,
                      Ipopt::Index&   nnz_jac_g,
                      Ipopt::Index&   nnz_h_lag,
                      IndexStyleEnum& index_style) override
    {
      n           = static_cast<Ipopt::Index>(model_.variableCount());
      m           = static_cast<Ipopt::Index>(model_.constraintCount());
      nnz_jac_g   = n * m;
      nnz_h_lag   = (n * (n + 1)) / 2;
      index_style = TNLP::C_STYLE;
      return true;
    }

    /// The latest negative model return; zero while none occurred.
    int fatalError() const
    {
      return fatal_;
    }

    bool get_bounds_info(Ipopt::Index,
                         Ipopt::Number* x_l,
                         Ipopt::Number* x_u,
                         Ipopt::Index,
                         Ipopt::Number* g_l,
                         Ipopt::Number* g_u) override
    {
      return accept(model_.variableBounds(x_l, x_u)) && accept(model_.constraintBounds(g_l, g_u));
    }

    bool get_starting_point(Ipopt::Index   n,
                            bool           init_x,
                            Ipopt::Number* x,
                            bool           init_z,
                            Ipopt::Number*,
                            Ipopt::Number*,
                            Ipopt::Index,
                            bool init_lambda,
                            Ipopt::Number*) override
    {
      if (init_z || init_lambda)
      {
        return false;
      }
      if (init_x)
      {
        for (Ipopt::Index i = 0; i < n; ++i)
        {
          x[i] = solution_[static_cast<size_t>(i)];
        }
      }
      return true;
    }

    bool eval_f(Ipopt::Index,
                const Ipopt::Number* x,
                bool,
                Ipopt::Number& obj_value) override
    {
      return accept(model_.evaluateObjective(x, obj_value));
    }

    bool eval_grad_f(Ipopt::Index,
                     const Ipopt::Number* x,
                     bool,
                     Ipopt::Number* grad_f) override
    {
      return accept(model_.evaluateGradient(x, grad_f));
    }

    bool eval_g(Ipopt::Index,
                const Ipopt::Number* x,
                bool,
                Ipopt::Index,
                Ipopt::Number* g) override
    {
      return accept(model_.evaluateConstraints(x, g));
    }

    bool eval_jac_g(Ipopt::Index         n,
                    const Ipopt::Number* x,
                    bool,
                    Ipopt::Index m,
                    Ipopt::Index,
                    Ipopt::Index*  iRow,
                    Ipopt::Index*  jCol,
                    Ipopt::Number* values) override
    {
      if (values == nullptr)
      {
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < m; ++row)
        {
          for (Ipopt::Index col = 0; col < n; ++col)
          {
            iRow[entry] = row;
            jCol[entry] = col;
            ++entry;
          }
        }
        return true;
      }
      return accept(model_.evaluateJacobian(x, values));
    }

    bool eval_h(Ipopt::Index         n,
                const Ipopt::Number* x,
                bool,
                Ipopt::Number obj_factor,
                Ipopt::Index,
                const Ipopt::Number* lambda,
                bool,
                Ipopt::Index,
                Ipopt::Index*  iRow,
                Ipopt::Index*  jCol,
                Ipopt::Number* values) override
    {
      if (!exact_hessian_)
      {
        return false;
      }

      if (values == nullptr)
      {
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < n; ++row)
        {
          for (Ipopt::Index col = 0; col <= row; ++col)
          {
            iRow[entry] = row;
            jCol[entry] = col;
            ++entry;
          }
        }
        return true;
      }

      if (!accept(model_.evaluateHessian(x, obj_factor, lambda, hessian_buffer_.data())))
      {
        return false;
      }

      const auto   columns = static_cast<size_t>(n);
      Ipopt::Index entry   = 0;
      for (Ipopt::Index row = 0; row < n; ++row)
      {
        for (Ipopt::Index col = 0; col <= row; ++col)
        {
          values[entry] = hessian_buffer_[static_cast<size_t>(row) * columns + static_cast<size_t>(col)];
          ++entry;
        }
      }
      return true;
    }

    void finalize_solution(Ipopt::SolverReturn,
                           Ipopt::Index         n,
                           const Ipopt::Number* x,
                           const Ipopt::Number*,
                           const Ipopt::Number*,
                           Ipopt::Index,
                           const Ipopt::Number*,
                           const Ipopt::Number*,
                           Ipopt::Number obj_value,
                           const Ipopt::IpoptData*,
                           Ipopt::IpoptCalculatedQuantities*) override
    {
      for (Ipopt::Index i = 0; i < n; ++i)
      {
        solution_[static_cast<size_t>(i)] = x[i];
      }
      stats_.final_objective = obj_value;
    }

    bool intermediate_callback(Ipopt::AlgorithmMode,
                               Ipopt::Index,
                               Ipopt::Number,
                               Ipopt::Number,
                               Ipopt::Number,
                               Ipopt::Number,
                               Ipopt::Number,
                               Ipopt::Number,
                               Ipopt::Number,
                               Ipopt::Number,
                               Ipopt::Index,
                               const Ipopt::IpoptData*,
                               Ipopt::IpoptCalculatedQuantities*) override
    {
      return fatal_ == 0;
    }

  private:
    /// Record a negative model return as fatal; true only on success.
    bool accept(int retval)
    {
      if (retval < 0)
      {
        fatal_ = retval;
      }
      return retval == 0;
    }

    const ModelT&       model_;
    std::vector<RealT>& solution_;
    StatsT&             stats_;
    bool                exact_hessian_;
    std::vector<RealT>  hessian_buffer_;
    int                 fatal_{0};
  };
} // namespace

namespace GridKit
{
  namespace Optimization
  {
    template <typename scalar_type, typename index_type>
    std::string IpoptSolver<scalar_type, index_type>::Stats::report() const
    {
      std::ostringstream stream;
      stream << "IpoptSolver: status " << status << ", iterations "
             << iterations << ", objective " << final_objective;
      return stream.str();
    }

    template <typename scalar_type, typename index_type>
    IpoptSolver<scalar_type, index_type>::IpoptSolver(
        OptimizationModel<ScalarT, IdxT>* model, Parameters params)
      : OptimizationSolver<ScalarT, IdxT>(model),
        params_(params)
    {
    }

    /**
     * @brief Solve the bound model.
     *
     * Normalizes @p x to a valid starting point (the model's default when
     * @p x is not sized to the problem), configures one Ipopt application,
     * and maps its return status onto the GridKit error convention:
     * 0 optimal, 1 acceptable, 2 iteration limit, negative otherwise. A
     * negative model-evaluator return aborts the iteration and is
     * reported as -3, honoring the OptimizationModel fatal contract.
     */
    template <typename scalar_type, typename index_type>
    int IpoptSolver<scalar_type, index_type>::solve(std::vector<RealT>& x)
    {
      const auto variables =
          static_cast<size_t>(this->model_->variableCount());
      if (x.size() != variables)
      {
        x.assign(variables, RealT{0});
        if (this->model_->startingPoint(x.data()) != 0)
        {
          return -1;
        }
      }
      stats_ = Stats{};

      // Honor the availability contract: a positive return from
      // evaluateHessian declares it unavailable, so the solve falls back
      // to the quasi-Newton approximation instead of failing later.
      bool exact_hessian = params_.hessian == HessianMode::EXACT;
      if (exact_hessian)
      {
        std::vector<RealT>       probe(variables * variables);
        const std::vector<RealT> multipliers(
            static_cast<size_t>(this->model_->constraintCount()));
        const int available = this->model_->evaluateHessian(
            x.data(), RealT{1}, multipliers.data(), probe.data());
        if (available < 0)
        {
          return -1;
        }
        exact_hessian = available == 0;
      }

      Ipopt::SmartPtr<Ipopt::IpoptApplication> application =
          IpoptApplicationFactory();
      application->Options()->SetNumericValue("tol", params_.tolerance);
      application->Options()->SetIntegerValue(
          "max_iter", static_cast<Ipopt::Index>(params_.max_iterations));
      application->Options()->SetIntegerValue(
          "print_level", static_cast<Ipopt::Index>(params_.print_level));
      application->Options()->SetStringValue("sb", "yes");
      // Variable bounds carry model contracts such as stability margins;
      // never let the interior point relax them.
      application->Options()->SetNumericValue("bound_relax_factor", 0.0);
      if (!exact_hessian)
      {
        application->Options()->SetStringValue("hessian_approximation",
                                               "limited-memory");
      }

      if (application->Initialize() != Ipopt::Solve_Succeeded)
      {
        return -1;
      }

      ModelAdapter<ScalarT, IdxT>* adapter =
          new ModelAdapter<ScalarT, IdxT>(*this->model_, x, stats_, exact_hessian);
      Ipopt::SmartPtr<Ipopt::TNLP> owner(adapter);

      const auto status = application->OptimizeTNLP(owner);
      stats_.status     = static_cast<int>(status);
      if (Ipopt::IsValid(application->Statistics()))
      {
        stats_.iterations = static_cast<IdxT>(
            application->Statistics()->IterationCount());
      }

      if (adapter->fatalError() != 0)
      {
        return -3;
      }

      switch (status)
      {
      case Ipopt::Solve_Succeeded:
        return 0;
      case Ipopt::Solved_To_Acceptable_Level:
        return 1;
      case Ipopt::Maximum_Iterations_Exceeded:
        return 2;
      default:
        return -2;
      }
    }

    template class IpoptSolver<double, long int>;
    template class IpoptSolver<double, size_t>;
  } // namespace Optimization
} // namespace GridKit

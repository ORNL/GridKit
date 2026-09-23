/**
 * @file OptimizationProblem.cpp
 * @brief Ipopt problem over a `Model::OptimizationEvaluator`.
 */

#include "OptimizationProblem.hpp"

#include <algorithm>

namespace AnalysisManager
{
  namespace IpoptInterface
  {
    template <typename scalar_type, typename index_type>
    OptimizationProblem<scalar_type, index_type>::OptimizationProblem(ModelT* model)
      : model_(model)
    {
    }

    template <typename scalar_type, typename index_type>
    bool OptimizationProblem<scalar_type, index_type>::get_nlp_info(Index&          n,
                                                                    Index&          m,
                                                                    Index&          nnz_jac_g,
                                                                    Index&          nnz_h_lag,
                                                                    IndexStyleEnum& index_style)
    {
      n           = static_cast<Index>(model_->size());
      m           = static_cast<Index>(model_->sizeConstraints());
      nnz_jac_g   = static_cast<Index>(model_->getCsrJacobian()->getNnz());
      nnz_h_lag   = static_cast<Index>(model_->getCsrHessian()->getNnz());
      index_style = C_STYLE;

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool OptimizationProblem<scalar_type, index_type>::get_bounds_info(Index   n,
                                                                       Number* x_l,
                                                                       Number* x_u,
                                                                       Index   m,
                                                                       Number* g_l,
                                                                       Number* g_u)
    {
      std::copy_n(model_->xLower().getData(), n, x_l);
      std::copy_n(model_->xUpper().getData(), n, x_u);
      std::copy_n(model_->gLower().getData(), m, g_l);
      std::copy_n(model_->gUpper().getData(), m, g_u);

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool OptimizationProblem<scalar_type, index_type>::get_starting_point(Index                    n,
                                                                          bool                     init_x,
                                                                          Number*                  x,
                                                                          bool                     init_z,
                                                                          [[maybe_unused]] Number* z_L,
                                                                          [[maybe_unused]] Number* z_U,
                                                                          [[maybe_unused]] Index   m,
                                                                          bool                     init_lambda,
                                                                          [[maybe_unused]] Number* lambda)
    {
      // Only the variables have a starting point
      if (!init_x || init_z || init_lambda)
      {
        return false;
      }

      std::copy_n(model_->x().getData(), n, x);

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool OptimizationProblem<scalar_type, index_type>::eval_f([[maybe_unused]] Index n,
                                                              const Number*          x,
                                                              [[maybe_unused]] bool  new_x,
                                                              Number&                obj_value)
    {
      setVariables(x);
      if (model_->evaluateObjective() != 0)
      {
        return false;
      }

      obj_value = model_->objective();

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool OptimizationProblem<scalar_type, index_type>::eval_grad_f(Index                 n,
                                                                   const Number*         x,
                                                                   [[maybe_unused]] bool new_x,
                                                                   Number*               grad_f)
    {
      setVariables(x);
      if (model_->evaluateGradient() != 0)
      {
        return false;
      }

      std::copy_n(model_->gradient().getData(), n, grad_f);

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool OptimizationProblem<scalar_type, index_type>::eval_g([[maybe_unused]] Index n,
                                                              const Number*          x,
                                                              [[maybe_unused]] bool  new_x,
                                                              Index                  m,
                                                              Number*                g)
    {
      setVariables(x);
      if (model_->evaluateConstraints() != 0)
      {
        return false;
      }

      std::copy_n(model_->g().getData(), m, g);

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool OptimizationProblem<scalar_type, index_type>::eval_jac_g([[maybe_unused]] Index n,
                                                                  const Number*          x,
                                                                  [[maybe_unused]] bool  new_x,
                                                                  [[maybe_unused]] Index m,
                                                                  Index                  nele_jac,
                                                                  Index*                 iRow,
                                                                  Index*                 jCol,
                                                                  Number*                values)
    {
      if (values == nullptr)
      {
        structure(*model_->getCsrJacobian(), iRow, jCol);
        return true;
      }

      setVariables(x);
      if (model_->evaluateJacobian() != 0)
      {
        return false;
      }

      std::copy_n(model_->getCsrJacobian()->getValues(), nele_jac, values);

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool OptimizationProblem<scalar_type, index_type>::eval_h([[maybe_unused]] Index n,
                                                              const Number*          x,
                                                              [[maybe_unused]] bool  new_x,
                                                              Number                 obj_factor,
                                                              [[maybe_unused]] Index m,
                                                              const Number*          lambda,
                                                              [[maybe_unused]] bool  new_lambda,
                                                              Index                  nele_hess,
                                                              Index*                 iRow,
                                                              Index*                 jCol,
                                                              Number*                values)
    {
      if (values == nullptr)
      {
        structure(*model_->getCsrHessian(), iRow, jCol);
        return true;
      }

      setVariables(x);
      if (model_->evaluateHessian(obj_factor, lambda) != 0)
      {
        return false;
      }

      std::copy_n(model_->getCsrHessian()->getValues(), nele_hess, values);

      return true;
    }

    /**
     * @brief Leave the final point in the model
     */
    template <typename scalar_type, typename index_type>
    void OptimizationProblem<scalar_type, index_type>::finalize_solution([[maybe_unused]] SolverReturn               status,
                                                                         [[maybe_unused]] Index                      n,
                                                                         const Number*                               x,
                                                                         [[maybe_unused]] const Number*              z_L,
                                                                         [[maybe_unused]] const Number*              z_U,
                                                                         [[maybe_unused]] Index                      m,
                                                                         [[maybe_unused]] const Number*              g,
                                                                         [[maybe_unused]] const Number*              lambda,
                                                                         [[maybe_unused]] Number                     obj_value,
                                                                         [[maybe_unused]] const IpoptData*           ip_data,
                                                                         [[maybe_unused]] IpoptCalculatedQuantities* ip_cq)
    {
      setVariables(x);
    }

    template <typename scalar_type, typename index_type>
    void OptimizationProblem<scalar_type, index_type>::setVariables(const Number* x)
    {
      std::copy_n(x, model_->size(), model_->x().getData());
      model_->x().setDataUpdated();
    }

    /**
     * @brief Coordinates of the CSR entries in storage order
     */
    template <typename scalar_type, typename index_type>
    void OptimizationProblem<scalar_type, index_type>::structure(CsrMatrixT& matrix, Index* iRow, Index* jCol)
    {
      const IdxT* row_ptrs = matrix.getRowData();
      const IdxT* cols     = matrix.getColData();
      for (IdxT row = 0; row < matrix.getNumRows(); ++row)
      {
        for (IdxT k = row_ptrs[row]; k < row_ptrs[row + 1]; ++k)
        {
          iRow[k] = static_cast<Index>(row);
          jCol[k] = static_cast<Index>(cols[k]);
        }
      }
    }

    // Available template instantiations
    template class OptimizationProblem<double, size_t>;
  } // namespace IpoptInterface
} // namespace AnalysisManager

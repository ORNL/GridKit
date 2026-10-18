/**
 * @file IpoptInterface.cpp
 * @brief Ipopt problem over a `Model::OptimizationEvaluator`.
 */

#include "IpoptInterface.hpp"

#include <algorithm>

namespace AnalysisManager
{
  namespace SteadyState
  {
    template <typename scalar_type, typename index_type>
    IpoptInterface<scalar_type, index_type>::IpoptInterface(ModelT* model)
      : model_(model)
    {
    }

    template <typename scalar_type, typename index_type>
    bool IpoptInterface<scalar_type, index_type>::get_nlp_info(ipindex&        n,
                                                               ipindex&        m,
                                                               ipindex&        nnz_jac_g,
                                                               ipindex&        nnz_h_lag,
                                                               IndexStyleEnum& index_style)
    {
      n           = static_cast<ipindex>(model_->size());
      m           = static_cast<ipindex>(model_->sizeEqualities() + model_->sizeInequalities());
      nnz_jac_g   = static_cast<ipindex>(model_->getCsrEqualityJacobian()->getNnz()
                                       + model_->getCsrInequalityJacobian()->getNnz());
      nnz_h_lag   = static_cast<ipindex>(model_->getCsrHessian()->getNnz());
      index_style = C_STYLE;

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool IpoptInterface<scalar_type, index_type>::get_bounds_info(ipindex   n,
                                                                  ipnumber* x_l,
                                                                  ipnumber* x_u,
                                                                  ipindex   m,
                                                                  ipnumber* g_l,
                                                                  ipnumber* g_u)
    {
      const ipindex m_c = static_cast<ipindex>(model_->sizeEqualities());

      std::copy_n(model_->xLower().getData(), n, x_l);
      std::copy_n(model_->xUpper().getData(), n, x_u);
      std::fill_n(g_l, m_c, 0.0);
      std::fill_n(g_u, m_c, 0.0);
      std::copy_n(model_->dLower().getData(), m - m_c, g_l + m_c);
      std::copy_n(model_->dUpper().getData(), m - m_c, g_u + m_c);

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool IpoptInterface<scalar_type, index_type>::get_starting_point(ipindex                    n,
                                                                     bool                       init_x,
                                                                     ipnumber*                  x,
                                                                     bool                       init_z,
                                                                     [[maybe_unused]] ipnumber* z_L,
                                                                     [[maybe_unused]] ipnumber* z_U,
                                                                     [[maybe_unused]] ipindex   m,
                                                                     bool                       init_lambda,
                                                                     [[maybe_unused]] ipnumber* lambda)
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
    bool IpoptInterface<scalar_type, index_type>::eval_f([[maybe_unused]] ipindex n,
                                                         const ipnumber*          x,
                                                         [[maybe_unused]] bool    new_x,
                                                         ipnumber&                obj_value)
    {
      setVariables(x);
      if (model_->evaluateObjective() != 0)
      {
        return false;
      }

      obj_value = model_->f();

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool IpoptInterface<scalar_type, index_type>::eval_grad_f(ipindex               n,
                                                              const ipnumber*       x,
                                                              [[maybe_unused]] bool new_x,
                                                              ipnumber*             grad_f)
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
    bool IpoptInterface<scalar_type, index_type>::eval_g([[maybe_unused]] ipindex n,
                                                         const ipnumber*          x,
                                                         [[maybe_unused]] bool    new_x,
                                                         ipindex                  m,
                                                         ipnumber*                g)
    {
      setVariables(x);
      if (model_->evaluateEqualities() != 0 || model_->evaluateInequalities() != 0)
      {
        return false;
      }

      const ipindex m_c = static_cast<ipindex>(model_->sizeEqualities());
      std::copy_n(model_->c().getData(), m_c, g);
      std::copy_n(model_->d().getData(), m - m_c, g + m_c);

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool IpoptInterface<scalar_type, index_type>::eval_jac_g([[maybe_unused]] ipindex n,
                                                             const ipnumber*          x,
                                                             [[maybe_unused]] bool    new_x,
                                                             [[maybe_unused]] ipindex m,
                                                             ipindex                  nele_jac,
                                                             ipindex*                 iRow,
                                                             ipindex*                 jCol,
                                                             ipnumber*                values)
    {
      if (values == nullptr)
      {
        const ipindex m_c   = static_cast<ipindex>(model_->sizeEqualities());
        const ipindex nnz_c = static_cast<ipindex>(model_->getCsrEqualityJacobian()->getNnz());
        structure(*model_->getCsrEqualityJacobian(), 0, iRow, jCol);
        structure(*model_->getCsrInequalityJacobian(), m_c, iRow + nnz_c, jCol + nnz_c);
        return true;
      }

      setVariables(x);
      if (model_->evaluateEqualityJacobian() != 0 || model_->evaluateInequalityJacobian() != 0)
      {
        return false;
      }

      const ipindex nnz_c = static_cast<ipindex>(model_->getCsrEqualityJacobian()->getNnz());
      std::copy_n(model_->getCsrEqualityJacobian()->getValues(), nnz_c, values);
      std::copy_n(model_->getCsrInequalityJacobian()->getValues(), nele_jac - nnz_c, values + nnz_c);

      return true;
    }

    template <typename scalar_type, typename index_type>
    bool IpoptInterface<scalar_type, index_type>::eval_h([[maybe_unused]] ipindex n,
                                                         const ipnumber*          x,
                                                         [[maybe_unused]] bool    new_x,
                                                         ipnumber                 obj_factor,
                                                         [[maybe_unused]] ipindex m,
                                                         const ipnumber*          lambda,
                                                         [[maybe_unused]] bool    new_lambda,
                                                         ipindex                  nele_hess,
                                                         ipindex*                 iRow,
                                                         ipindex*                 jCol,
                                                         ipnumber*                values)
    {
      if (values == nullptr)
      {
        structure(*model_->getCsrHessian(), 0, iRow, jCol);
        return true;
      }

      setVariables(x);
      const ipindex m_c = static_cast<ipindex>(model_->sizeEqualities());
      if (model_->evaluateHessian(obj_factor, lambda, lambda + m_c) != 0)
      {
        return false;
      }

      std::copy_n(model_->getCsrHessian()->getValues(), nele_hess, values);

      return true;
    }

    template <typename scalar_type, typename index_type>
    void IpoptInterface<scalar_type, index_type>::finalize_solution([[maybe_unused]] SolverReturn               status,
                                                                    [[maybe_unused]] ipindex                    n,
                                                                    const ipnumber*                             x,
                                                                    [[maybe_unused]] const ipnumber*            z_L,
                                                                    [[maybe_unused]] const ipnumber*            z_U,
                                                                    [[maybe_unused]] ipindex                    m,
                                                                    [[maybe_unused]] const ipnumber*            g,
                                                                    [[maybe_unused]] const ipnumber*            lambda,
                                                                    [[maybe_unused]] ipnumber                   obj_value,
                                                                    [[maybe_unused]] const IpoptData*           ip_data,
                                                                    [[maybe_unused]] IpoptCalculatedQuantities* ip_cq)
    {
      setVariables(x);
      model_->evaluateObjective();
    }

    template <typename scalar_type, typename index_type>
    void IpoptInterface<scalar_type, index_type>::setVariables(const ipnumber* x)
    {
      std::copy_n(x, model_->size(), model_->x().getData());
      model_->x().setDataUpdated();
    }

    template <typename scalar_type, typename index_type>
    void IpoptInterface<scalar_type, index_type>::structure(CsrMatrixT& matrix,
                                                            ipindex     row_offset,
                                                            ipindex*    iRow,
                                                            ipindex*    jCol)
    {
      const IdxT* row_ptrs = matrix.getRowData();
      const IdxT* cols     = matrix.getColData();
      for (IdxT row = 0; row < matrix.getNumRows(); ++row)
      {
        for (IdxT k = row_ptrs[row]; k < row_ptrs[row + 1]; ++k)
        {
          iRow[k] = row_offset + static_cast<ipindex>(row);
          jCol[k] = static_cast<ipindex>(cols[k]);
        }
      }
    }

    // Available template instantiations
    template class IpoptInterface<double, size_t>;
  } // namespace SteadyState
} // namespace AnalysisManager

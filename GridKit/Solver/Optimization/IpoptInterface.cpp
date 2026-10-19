/**
 * @file IpoptInterface.cpp
 * @brief Ipopt problem over a `Model::OptimizationEvaluator`.
 */

#include "IpoptInterface.hpp"

#include <algorithm>

#include <GridKit/Utilities/Logger/Logger.hpp>

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
      if (model_->getCsrJacobian() == nullptr || model_->getCsrHessian() == nullptr)
      {
        GridKit::Utilities::Logger::error() << "IpoptInterface: allocate() the model before solving it\n";
        return false;
      }

      n           = static_cast<ipindex>(model_->size());
      m           = static_cast<ipindex>(model_->sizeConstraints());
      nnz_jac_g   = static_cast<ipindex>(model_->getCsrJacobian()->getNnz());
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
      std::copy_n(model_->xLower().getData(), n, x_l);
      std::copy_n(model_->xUpper().getData(), n, x_u);
      std::copy_n(model_->gLower().getData(), m, g_l);
      std::copy_n(model_->gUpper().getData(), m, g_u);

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
      // Only the variables have initial values
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
      copyInto(model_->x(), x);
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
      copyInto(model_->x(), x);
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
      copyInto(model_->x(), x);
      if (model_->evaluateConstraints() != 0)
      {
        return false;
      }

      std::copy_n(model_->g().getData(), m, g);

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
        structure(*model_->getCsrJacobian(), iRow, jCol);
        return true;
      }

      copyInto(model_->x(), x);
      if (model_->evaluateJacobian() != 0)
      {
        return false;
      }

      std::copy_n(model_->getCsrJacobian()->getValues(), nele_jac, values);

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
        structure(*model_->getCsrHessian(), iRow, jCol);
        return true;
      }

      copyInto(model_->x(), x);
      copyInto(model_->lambda(), lambda);
      if (model_->evaluateHessian(obj_factor) != 0)
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
                                                                    const ipnumber*                             z_L,
                                                                    const ipnumber*                             z_U,
                                                                    [[maybe_unused]] ipindex                    m,
                                                                    [[maybe_unused]] const ipnumber*            g,
                                                                    const ipnumber*                             lambda,
                                                                    [[maybe_unused]] ipnumber                   obj_value,
                                                                    [[maybe_unused]] const IpoptData*           ip_data,
                                                                    [[maybe_unused]] IpoptCalculatedQuantities* ip_cq)
    {
      copyInto(model_->x(), x);
      copyInto(model_->lambda(), lambda);
      copyInto(model_->zLower(), z_L);
      copyInto(model_->zUpper(), z_U);
      model_->evaluateObjective();
    }

    template <typename scalar_type, typename index_type>
    template <typename VectorT>
    void IpoptInterface<scalar_type, index_type>::copyInto(VectorT& vector, const ipnumber* values)
    {
      std::copy_n(values, vector.getSize(), vector.getData());
      vector.setDataUpdated();
    }

    template <typename scalar_type, typename index_type>
    void IpoptInterface<scalar_type, index_type>::structure(CsrMatrixT& matrix, ipindex* iRow, ipindex* jCol)
    {
      const IdxT* row_ptrs = matrix.getRowData();
      const IdxT* cols     = matrix.getColData();
      for (IdxT row = 0; row < matrix.getNumRows(); ++row)
      {
        for (IdxT k = row_ptrs[row]; k < row_ptrs[row + 1]; ++k)
        {
          iRow[k] = static_cast<ipindex>(row);
          jCol[k] = static_cast<ipindex>(cols[k]);
        }
      }
    }

    // Available template instantiations
    template class IpoptInterface<ipnumber, size_t>;
  } // namespace SteadyState
} // namespace AnalysisManager

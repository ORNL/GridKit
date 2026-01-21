/**
 * @file DhDy.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/ModelWrappers.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/COO_Matrix.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: dh/dy
       *
       * @tparam ModelT - model type
       * @tparam MemberFunctions - member function parameter key
       * @tparam ScalarT - scalar data type
       * @tparam IdxT - matrix index data type
       */
      template <typename ModelT, MemberFunctions function, class ScalarT, typename IdxT>
      struct DhDy
      {
        using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;
        using MatrixT = GridKit::LinearAlgebra::COO_Matrix<RealT, IdxT>;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_res - Number of residual functions
         * @param[in] n_var - Number of independent variables
         * @param[in] res_indices - Global residual indices
         * @param[in] var_indices - Global variable indices
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[in,out] jac - Jacobian
         */
        static void eval(ModelT*                  model,
                         size_t                   n_res,
                         size_t                   n_var,
                         const std::vector<IdxT>& res_indices,
                         const std::vector<IdxT>& var_indices,
                         ScalarT*                 y,
                         ScalarT*                 yp,
                         ScalarT*                 wb,
                         MatrixT&                 jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            std::vector<Triple<ScalarT>> triplets; //< @todo: Update once sparse storage format changes
            std::vector<ScalarT>         elementary_v(n_var);
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, var_i, &triplets);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Core automatic differentiaation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function, ScalarT>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_dup,
                                     y,
                                     output,
                                     enzyme_const,
                                     yp,
                                     enzyme_const,
                                     wb,
                                     enzyme_dupnoneed,
                                     elementary_v.data(),
                                     d_output);
            }

            // Store result
            std::vector<IdxT>    ctemp{};
            std::vector<IdxT>    rtemp{};
            std::vector<ScalarT> valtemp{};
            for (auto& tup : triplets)
            {
              rtemp.push_back(res_indices[static_cast<size_t>(tup.row)]);
              ctemp.push_back(var_indices[static_cast<size_t>(tup.col)]);
              valtemp.push_back(tup.val);
            }
            jac.setValues(rtemp, ctemp, valtemp); //< @todo: Update once sparse storage format changes
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit

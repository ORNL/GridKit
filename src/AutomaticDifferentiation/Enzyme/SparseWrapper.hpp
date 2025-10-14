/**
 * @file SparseWrapper.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <vector>

/**
 * @brief Enzyme constants for activity analysis
 *
 */
extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_dupnoneed;

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Model member function parameter keys
       *
       */
      enum class MemberFunctions
      {
        InternalResidual,
        InternalResidualWithExternal,
        BusResidual,
        BusResidual11, //< Special case for branches that are connected to two buses
        BusResidual12, //< Special case for branches that are connected to two buses
        BusResidual21, //< Special case for branches that are connected to two buses
        BusResidual22  //< Special case for branches that are connected to two buses
      };

      /**
       * @brief Residual wrapper around residual methods inside model classes
       *
       * @tparam ModelT - model type
       * @tparam MemberFunctions - member function parameter key
       * @tparam ScalarT - scalar data type
       *
       * @param[in] model - Pointer to the model to be differentiated
       * @param[in] y - Internal variables
       * @param[in] yp - Internal variable derivatives
       * @param[in] wb - Coupling (bus) variables
       * @param[in] we - External variables
       * @param[out] f - Internal residual
       * @param[out] h - Coupling (bus) residual
       */
      template <typename ModelT, MemberFunctions function, typename ScalarT>
      struct Wrapper
      {
        static void eval(ModelT* model, ScalarT* y, ScalarT* yp, ScalarT* wb, ScalarT* f)
        {
          model->evaluateInternalResidual(y, yp, wb, f);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for InternalResidualWithExternal
       *
       */
      template <typename ModelT, typename ScalarT>
      struct Wrapper<ModelT, MemberFunctions::InternalResidualWithExternal, ScalarT>
      {
        static void eval(ModelT* model, ScalarT* y, ScalarT* yp, ScalarT* wb, ScalarT* we, ScalarT* h)
        {
          model->evaluateInternalResidualWithExternal(y, yp, wb, we, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual
       *
       */
      template <typename ModelT, typename ScalarT>
      struct Wrapper<ModelT, MemberFunctions::BusResidual, ScalarT>
      {
        static void eval(ModelT* model, ScalarT* y, ScalarT* yp, ScalarT* wb, ScalarT* h)
        {
          model->evaluateBusResidual(y, yp, wb, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual11
       *
       */
      template <typename ModelT, typename ScalarT>
      struct Wrapper<ModelT, MemberFunctions::BusResidual11, ScalarT>
      {
        static void eval(ModelT* model, ScalarT* y, ScalarT* yp, ScalarT* wb, ScalarT* h)
        {
          model->evaluateBusResidual11(y, yp, wb, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual12
       *
       */
      template <typename ModelT, typename ScalarT>
      struct Wrapper<ModelT, MemberFunctions::BusResidual12, ScalarT>
      {
        static void eval(ModelT* model, ScalarT* y, ScalarT* yp, ScalarT* wb, ScalarT* h)
        {
          model->evaluateBusResidual12(y, yp, wb, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual21
       *
       */
      template <typename ModelT, typename ScalarT>
      struct Wrapper<ModelT, MemberFunctions::BusResidual21, ScalarT>
      {
        static void eval(ModelT* model, ScalarT* y, ScalarT* yp, ScalarT* wb, ScalarT* h)
        {
          model->evaluateBusResidual21(y, yp, wb, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual22
       *
       */
      template <typename ModelT, typename ScalarT>
      struct Wrapper<ModelT, MemberFunctions::BusResidual22, ScalarT>
      {
        static void eval(ModelT* model, ScalarT* y, ScalarT* yp, ScalarT* wb, ScalarT* h)
        {
          model->evaluateBusResidual22(y, yp, wb, h);
        }
      };

      /**
       * @brief Enzyme fwddiff template
       *
       * @tparam T - return type
       * @tparam ModelT - model type
       */
      template <typename T, typename... ModelT>
      extern T __enzyme_fwddiff(void*, ModelT...) noexcept;

      /**
       * @brief Enzyme todense template
       *
       * @tparam T - return type
       */
      template <typename T>
      extern T __enzyme_todense(void*...) noexcept;

      /**
       * @brief Enzyme sparse storage in triplet format
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      struct Triple
      {
        size_t  row;
        size_t  col;
        ScalarT val;
        Triple(Triple&&) = default;

        Triple(size_t row, size_t col, ScalarT val)
          : row(row),
            col(col),
            val(val)
        {
        }
      };

      /**
       * @brief Enzyme sparse accumulation for float
       *
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_storeflt(size_t row, size_t col, float val, std::vector<Triple<float>>& triplets)
      {
        triplets.emplace_back(row, col, val);
      }

      /**
       * @brief Enzyme sparse accumulation for double
       *
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_storedbl(size_t row, size_t col, double val, std::vector<Triple<double>>& triplets)
      {
        triplets.emplace_back(row, col, val);
      }

      /**
       * @brief Enzyme sparse store
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      __attribute__((always_inline)) static void sparse_store(ScalarT val, size_t idx, size_t i, std::vector<Triple<ScalarT>>& triplets)
      {
        if (val == 0.0)
          return;
        idx /= sizeof(ScalarT);
        if constexpr (sizeof(ScalarT) == 4)
          inner_storeflt(idx, i, val, triplets);
        else
          inner_storedbl(idx, i, val, triplets);
      }

      /**
       * @brief Enzyme sparse load
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      __attribute__((always_inline)) static ScalarT sparse_load(size_t, size_t, std::vector<Triple<ScalarT>>&)
      {
        return 0.0;
      }

      /**
       * @brief Enzyme identity store
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      __attribute__((always_inline)) static void ident_store(ScalarT, size_t, size_t)
      {
        assert(0 && "should never store");
      }

      /**
       * @brief Enzyme identity load
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      __attribute__((always_inline)) static ScalarT ident_load(size_t idx, size_t i)
      {
        idx /= sizeof(ScalarT);
        return (ScalarT) (idx == i);
      }

      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: Internal Jacobian
       *
       * @tparam ModelT - model type
       * @tparam MemberFunctions - member function parameter key
       * @tparam ScalarT - scalar data type
       * @tparam IdxT - matrix index data type
       *
       * @param[in] model - Pointer to the model to be differentiated
       * @param[in] n_res - Number of residual functions
       * @param[in] n_var - Number of independent variables
       * @param[in] res_indices - Map from local residual indices to global indices
       * @param[in] var_indices - Map from local variable indices to global indices
       * @param[in] y - Internal variables
       * @param[in] yp - Internal variable derivatives
       * @param[in] wb - Coupling (bus) variables
       * @param[in] we - External variables
       * @param[in,out] jac - Jacobian
       */
      template <typename ModelT, MemberFunctions function, class ScalarT, typename IdxT>
      struct InternalJacobian
      {
        static void eval(ModelT*                                            model,
                         size_t                                             n_res,
                         size_t                                             n_var,
                         const std::map<IdxT, IdxT>&                        res_indices,
                         const std::map<IdxT, IdxT>&                        var_indices,
                         ScalarT*                                           y,
                         ScalarT*                                           yp,
                         ScalarT*                                           wb,
                         GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            std::vector<Triple<ScalarT>> triplets;
            std::vector<ScalarT>         elementary_v(n_res);
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, var_i, &triplets);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Autodiff
              __enzyme_fwddiff<void>((void*) Wrapper<ModelT, function, ScalarT>::eval,
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
              rtemp.push_back(res_indices.at(static_cast<IdxT>(tup.row)));
              ctemp.push_back(var_indices.at(static_cast<IdxT>(tup.col)));
              valtemp.push_back(tup.val);
            }
            jac.setValues(rtemp, ctemp, valtemp); //< @todo: Update once sparse storage format changes
          }
        }
      };

      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: Internal Jacobian with external variables
       *
       */
      template <typename ModelT, MemberFunctions function, class ScalarT, typename IdxT>
      struct InternalJacobianWithExternal
      {
        static void eval(ModelT*                                            model,
                         size_t                                             n_res,
                         size_t                                             n_var,
                         const std::map<IdxT, IdxT>&                        res_indices,
                         const std::map<IdxT, IdxT>&                        var_indices,
                         ScalarT*                                           y,
                         ScalarT*                                           yp,
                         ScalarT*                                           wb,
                         ScalarT*                                           we,
                         GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            std::vector<Triple<ScalarT>> triplets;
            std::vector<ScalarT>         elementary_v(n_res);
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, var_i, &triplets);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Autodiff
              __enzyme_fwddiff<void>((void*) Wrapper<ModelT, function, ScalarT>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_dup,
                                     y,
                                     output,
                                     enzyme_const,
                                     yp,
                                     enzyme_const,
                                     wb,
                                     enzyme_const,
                                     we,
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
              rtemp.push_back(res_indices.at(static_cast<IdxT>(tup.row)));
              ctemp.push_back(var_indices.at(static_cast<IdxT>(tup.col)));
              valtemp.push_back(tup.val);
            }
            jac.setValues(rtemp, ctemp, valtemp); //< @todo: Update once sparse storage format changes
          }
        }
      };

      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: External Jacobian
       *
       */
      template <typename ModelT, MemberFunctions function, class ScalarT, typename IdxT>
      struct ExternalJacobian
      {
        static void eval(ModelT*                                            model,
                         size_t                                             n_res,
                         size_t                                             n_var,
                         const std::map<IdxT, IdxT>&                        res_indices,
                         const std::map<IdxT, IdxT>&                        var_indices,
                         ScalarT*                                           y,
                         ScalarT*                                           yp,
                         ScalarT*                                           wb,
                         ScalarT*                                           we,
                         GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            std::vector<Triple<ScalarT>> triplets;
            std::vector<ScalarT>         elementary_v(n_res);
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, var_i, &triplets);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Autodiff
              __enzyme_fwddiff<void>((void*) Wrapper<ModelT, function, ScalarT>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_const,
                                     y,
                                     enzyme_const,
                                     yp,
                                     enzyme_const,
                                     wb,
                                     enzyme_dup,
                                     we,
                                     output,
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
              IdxT row = res_indices.at(static_cast<IdxT>(tup.row));
              IdxT col = var_indices.at(static_cast<IdxT>(tup.col));
              if (col != static_cast<IdxT>(-1))
              {
                rtemp.push_back(row);
                ctemp.push_back(col);
                valtemp.push_back(tup.val);
              }
            }
            jac.setValues(rtemp, ctemp, valtemp); //< @todo: Update once sparse storage format changes
          }
        }
      };

      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: Bus Jacobian
       *
       */
      template <typename ModelT, MemberFunctions function, class ScalarT, typename IdxT>
      struct BusJacobian
      {
        static void eval(ModelT*                                            model,
                         size_t                                             n_res,
                         size_t                                             n_var,
                         const std::map<IdxT, IdxT>&                        res_indices,
                         const std::map<IdxT, IdxT>&                        var_indices,
                         ScalarT*                                           y,
                         ScalarT*                                           yp,
                         ScalarT*                                           wb,
                         GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            std::vector<Triple<ScalarT>> triplets;
            std::vector<ScalarT>         elementary_v(n_res);
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, var_i, &triplets);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Autodiff
              __enzyme_fwddiff<void>((void*) Wrapper<ModelT, function, ScalarT>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_const,
                                     y,
                                     enzyme_const,
                                     yp,
                                     enzyme_dup,
                                     wb,
                                     output,
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
              rtemp.push_back(res_indices.at(static_cast<IdxT>(tup.row)));
              ctemp.push_back(var_indices.at(static_cast<IdxT>(tup.col)));
              valtemp.push_back(tup.val);
            }
            jac.axpy(1.0, rtemp, ctemp, valtemp); //< @todo: Update once sparse storage format changes
          }
        }
      };

      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: Branch Jacobian
       *
       */
      template <typename ModelT, MemberFunctions function, class ScalarT, typename IdxT>
      struct BranchJacobian
      {
        static void eval(ModelT*                                            model,
                         size_t                                             n_res,
                         size_t                                             n_var,
                         const std::map<IdxT, IdxT>&                        res_indices,
                         const std::map<IdxT, IdxT>&                        var_indices,
                         ScalarT*                                           y,
                         ScalarT*                                           yp,
                         ScalarT*                                           wb,
                         GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            std::vector<Triple<ScalarT>> triplets;
            std::vector<ScalarT>         elementary_v(n_res);
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, var_i, &triplets);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Autodiff
              __enzyme_fwddiff<void>((void*) Wrapper<ModelT, function, ScalarT>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_const,
                                     y,
                                     enzyme_const,
                                     yp,
                                     enzyme_dup,
                                     wb,
                                     output,
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
              rtemp.push_back(res_indices.at(static_cast<IdxT>(tup.row)));
              ctemp.push_back(var_indices.at(static_cast<IdxT>(tup.col)));
              valtemp.push_back(tup.val);
            }
            jac.setValues(rtemp, ctemp, valtemp); //< @todo: Update once sparse storage format changes
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit

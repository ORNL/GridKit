/**
 * @file ModelWrappers.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <vector>

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
        ExternalResidual
      };

      /**
       * @brief Template definition for wrapper around residual methods inside model classes
       *
       * @tparam ModelT - model type
       * @tparam MemberFunctions - member function parameter key
       *
       */
      template <typename ModelT, MemberFunctions function>
      struct ModelWrapper
      {
      };

      /**
       * @brief Residual wrapper partial template specialization for InternalResidual
       *
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::InternalResidual>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] y_ext - External variables
         * @param[out] f - Internal residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* y_ext,
                         ScalarT*       f)
        {
          model->evaluateInternalResidual(y, yp, y_ext, f);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for ExternalResidual
       *
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::ExternalResidual>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] y_ext - External variables
         * @param[out] f_ext - External residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* y_ext,
                         ScalarT*       f_ext)
        {
          model->evaluateExternalResidual(y, yp, y_ext, f_ext);
        }
      };

    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit

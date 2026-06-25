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
        InternalResidualWithSignal,
        BusResidual,
        BusResidual11, //< Special case for branches that are connected to two buses
        BusResidual12, //< Special case for branches that are connected to two buses
        BusResidual21, //< Special case for branches that are connected to two buses
        BusResidual22  //< Special case for branches that are connected to two buses
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
         * @param[in] wb - Bus variables
         * @param[out] f - Internal residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* wb,
                         ScalarT*       f)
        {
          model->evaluateInternalResidual(y, yp, wb, f);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for InternalResidualWithSignal
       *
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::InternalResidualWithSignal>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[in] ws - Signal variables
         * @param[out] f - Internal residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* wb,
                         const ScalarT* ws,
                         ScalarT*       f)
        {
          model->evaluateInternalResidual(y, yp, wb, ws, f);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual
       *
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::BusResidual>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[out] h - Bus residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* wb,
                         ScalarT*       h)
        {
          model->evaluateBusResidual(y, yp, wb, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual11 (branch member function)
       *
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::BusResidual11>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[out] h - Bus residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* wb,
                         ScalarT*       h)
        {
          model->evaluateBusResidual11(y, yp, wb, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual12 (branch member function)
       *
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::BusResidual12>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[out] h - Bus residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* wb,
                         ScalarT*       h)
        {
          model->evaluateBusResidual12(y, yp, wb, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual21 (branch member function)
       *
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::BusResidual21>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function>::eval,
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[out] h - Bus residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* wb,
                         ScalarT*       h)
        {
          model->evaluateBusResidual21(y, yp, wb, h);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for BusResidual22 (branch member function)
       *
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::BusResidual22>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[out] h - Bus residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* wb,
                         ScalarT*       h)
        {
          model->evaluateBusResidual22(y, yp, wb, h);
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit

/**
 * @file ModelWrappers.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <cstddef>
#include <vector>

#include <GridKit/Model/EvaluatorTraits.hpp>
#include <GridKit/Utilities/Enum.hpp>

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
        BusResidual22, //< Special case for branches that are connected to two buses
        Objective,     //< Optimization objective
        Constraints    //< Optimization internal constraints and contributions at the ports
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
        using ScalarT = typename GridKit::Model::EvaluatorTraits<ModelT>::ScalarT;

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
        using ScalarT = typename GridKit::Model::EvaluatorTraits<ModelT>::ScalarT;

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
        using ScalarT = typename GridKit::Model::EvaluatorTraits<ModelT>::ScalarT;

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
        using ScalarT = typename GridKit::Model::EvaluatorTraits<ModelT>::ScalarT;

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
        using ScalarT = typename GridKit::Model::EvaluatorTraits<ModelT>::ScalarT;

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
        using ScalarT = typename GridKit::Model::EvaluatorTraits<ModelT>::ScalarT;

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
        using ScalarT = typename GridKit::Model::EvaluatorTraits<ModelT>::ScalarT;

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

      /**
       * @brief Wrapper partial template specialization for an optimization objective
       *
       * Local variables are the internal variables, then the external ones.
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::Objective>
      {
        /// Number of internal variables
        static constexpr std::size_t NUM_INTERNAL_INPUTS = Utilities::enum_size<typename ModelT::InternalVariablesT>();

        /// Number of local variables
        static constexpr std::size_t NUM_INPUTS = NUM_INTERNAL_INPUTS + Utilities::enum_size<typename ModelT::ExternalVariablesT>();

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] x - Local variables
         * @return Objective contribution
         */
        template <typename ScalarT>
        [[gnu::always_inline]] inline static ScalarT eval(const ModelT* model, const ScalarT* x)
        {
          return model->objective(x, x + NUM_INTERNAL_INPUTS);
        }
      };

      /**
       * @brief Wrapper partial template specialization for optimization constraints
       *
       * Local contributions are those to the internal constraints, then those at the ports.
       */
      template <typename ModelT>
      struct ModelWrapper<ModelT, MemberFunctions::Constraints>
      {
        /// Number of internal variables
        static constexpr std::size_t NUM_INTERNAL_INPUTS = Utilities::enum_size<typename ModelT::InternalVariablesT>();

        /// Number of internal constraints
        static constexpr std::size_t NUM_INTERNAL_OUTPUTS = Utilities::enum_size<typename ModelT::InternalConstraintsT>();

        /// Number of local constraint contributions
        static constexpr std::size_t NUM_OUTPUTS = NUM_INTERNAL_OUTPUTS + Utilities::enum_size<typename ModelT::ExternalConstraintsT>();

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] x - Local variables
         * @param[out] g - Local constraint contributions
         */
        template <typename ScalarT>
        [[gnu::always_inline]] inline static void eval(const ModelT* model, const ScalarT* x, ScalarT* g)
        {
          model->internalConstraints(x, x + NUM_INTERNAL_INPUTS, g);
          model->externalConstraints(x, x + NUM_INTERNAL_INPUTS, g + NUM_INTERNAL_OUTPUTS);
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit

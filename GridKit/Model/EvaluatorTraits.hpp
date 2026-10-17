/**
 * @file EvaluatorTraits.hpp
 * @brief Scalar and index types of a model, taken from its Evaluator base.
 */

#pragma once

#include <utility>

#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief Exposes the scalar, index, and real types of any model that
     * publicly derives from Evaluator.
     *
     * Model classes keep their own type aliases private, so generic code
     * (for example the Enzyme Jacobian evaluators) obtains the types through
     * this trait instead of through the model class.
     *
     * @tparam model_type Model type derived from Evaluator
     */
    template <typename model_type>
    struct EvaluatorTraits
    {
    private:
      template <typename scalar_type, typename index_type>
      static Evaluator<scalar_type, index_type> probe(const Evaluator<scalar_type, index_type>&);

    public:
      using EvaluatorT = decltype(probe(std::declval<const model_type&>()));
      using ScalarT    = typename EvaluatorT::ScalarT;
      using IdxT       = typename EvaluatorT::IdxT;
      using RealT      = typename EvaluatorT::RealT;
    };
  } // namespace Model
} // namespace GridKit

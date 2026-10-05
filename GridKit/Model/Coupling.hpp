#pragma once

#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  namespace Model
  {
    /// A model input, held at value or ramped linearly from start.
    template <class ScalarT>
    struct Input
    {
      using RealT = typename ScalarTraits<ScalarT>::RealT;

      ScalarT value{};
      ScalarT rate{};
      RealT   start{};

      ScalarT at(RealT t) const
      {
        return value + rate * (t - start);
      }
    };

    /**
     * @brief An input of one model set from an entry of another model's state.
     *
     * A partitioned solver sets the input from state(source)[index] before the
     * coupled model is advanced.
     */
    template <class ScalarT, typename IdxT>
    struct Coupling
    {
      const Evaluator<ScalarT, IdxT>* source; ///< Model whose state holds the value
      IdxT                            index;  ///< Index in the source state
      Input<ScalarT>*                 input;  ///< Input of the coupled model
    };
  } // namespace Model
} // namespace GridKit

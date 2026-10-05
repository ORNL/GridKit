#pragma once

#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief An input of one model set from an entry of another model's state.
     *
     * A partitioned solver copies state(source)[index] into *value before the
     * coupled model is advanced.
     */
    template <class ScalarT, typename IdxT>
    struct Coupling
    {
      const Evaluator<ScalarT, IdxT>* source; ///< Model whose state holds the value
      IdxT                            index;  ///< Index in the source state
      ScalarT*                        value;  ///< Input of the coupled model
    };
  } // namespace Model
} // namespace GridKit

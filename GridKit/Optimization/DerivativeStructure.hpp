#pragma once

#include <cstddef>
#include <span>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Structurally present entry in a local constraint Jacobian.
     *
     * Component descriptors are ordered by variable and then constraint. The
     * matching derivative value always occupies the same descriptor position.
     * Local-to-global variable maps must contain unique indices and act only as
     * gathers. Constraint maps are additive scatters. Signs and scaling belong
     * in the primal equations, not in either map. Descriptors are unique and
     * must include every entry that can be nonzero at any admissible state.
     */
    template <typename index_type>
    struct LocalJacobianEntry
    {
      index_type constraint;
      index_type variable;
    };

    /**
     * @brief Structurally present entry in a local Lagrangian Hessian.
     *
     * Component descriptors contain the lower triangle and are ordered by
     * column and then row. They are unique and include every entry that can be
     * nonzero at any admissible state. Numerical zeros remain structurally
     * present.
     */
    template <typename index_type>
    struct LocalHessianEntry
    {
      index_type row;
      index_type column;
    };

    /**
     * @brief Matrix coordinate after local indices have been mapped globally.
     */
    template <typename index_type>
    struct SparseEntry
    {
      index_type row;
      index_type column;
    };

    /**
     * @brief Return the canonical lower-triangular coordinate for an entry.
     */
    template <typename index_type>
    constexpr SparseEntry<index_type> lowerTriangle(index_type row, index_type column)
    {
      if (row < column)
      {
        return {column, row};
      }
      return {row, column};
    }

    /**
     * @brief Check the uniqueness required of a local-to-global variable map.
     */
    template <typename index_type>
    constexpr bool hasUniqueIndices(std::span<const index_type> indices)
    {
      for (std::size_t first = 0; first < indices.size(); ++first)
      {
        for (std::size_t second = first + 1; second < indices.size(); ++second)
        {
          if (indices[first] == indices[second])
          {
            return false;
          }
        }
      }
      return true;
    }
  } // namespace Optimization
} // namespace GridKit

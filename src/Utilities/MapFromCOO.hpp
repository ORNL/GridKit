
/**
 * @file MapToCOO.hpp
 *
 * @author Nicholson Koukpaizan <koukpaizannk@ornl.gov>, ORNL
 * @todo This should go away once we settle on a sparse format
 *
 */
#pragma once

#include <AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <LinearAlgebra/SparseMatrix/COO_Matrix.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    std::vector<DependencyTracking::Variable::DependencyMap> MapFromCOO(LinearAlgebra::COO_Matrix<ScalarT, IdxT> matrix)
    {
      std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> matrix_entries = matrix.getEntries();
      const auto [rows, columns, values]                                                       = matrix_entries;

      std::tuple<IdxT, IdxT> matrix_dimensions = matrix.getDimensions();
      const auto [n_rows, n_columns]           = matrix_dimensions;

      std::vector<DependencyTracking::Variable::DependencyMap> dependencies(n_rows);

      for (IdxT i = 0; i < rows.size(); ++i)
      {
        dependencies[rows[i]].insert(std::make_pair(columns[i], values[i]));
      }

      return dependencies;
    }
  } // namespace Testing
} // namespace GridKit

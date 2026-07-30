
/**
 * @file MapFromCsr.hpp
 *
 * @author Nicholson Koukpaizan <koukpaizannk@ornl.gov>, ORNL
 * @todo Move to a more permanent location
 *
 */
#pragma once

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename real_type, typename index_type>
    std::vector<DependencyTracking::Variable::DependencyMap> MapFromCsr(LinearAlgebra::CsrMatrix<real_type, index_type>* matrix)
    {
      using RealT = real_type;
      using IdxT  = index_type;

      IdxT*  row_data    = matrix->getRowData();
      IdxT*  column_data = matrix->getColData();
      RealT* values      = matrix->getValues();
      IdxT   n_rows      = matrix->getNumRows();

      std::vector<DependencyTracking::Variable::DependencyMap> dependencies(n_rows);

      for (IdxT i = 0; i < n_rows; ++i)
      {
        for (IdxT j = row_data[i]; j < row_data[i + 1]; ++j)
        {
          dependencies[i].insert(std::make_pair(column_data[j], values[j]));
        }
      }

      return dependencies;
    }
  } // namespace Testing
} // namespace GridKit

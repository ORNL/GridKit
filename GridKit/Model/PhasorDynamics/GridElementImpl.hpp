#pragma once

#include <GridKit/Model/PhasorDynamics/GridElement.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename ScalarP, typename IdxP>
    GridElement<ScalarP, IdxP>::~GridElement()
    {
      if (J_rows_buffer_ != nullptr)
      {
        delete[] J_rows_buffer_;
        delete[] J_cols_buffer_;
        delete[] J_vals_buffer_;
        J_rows_buffer_ = nullptr;
        J_cols_buffer_ = nullptr;
        J_vals_buffer_ = nullptr;
      }
    }

    template <typename ScalarP, typename IdxP>
    void GridElement<ScalarP, IdxP>::setTolerances(RealT& rtol, RealT& atol) const
    {
      rtol = rel_tol_;
      atol = abs_tol_;
    }

    template <typename ScalarP, typename IdxP>
    void GridElement<ScalarP, IdxP>::setMaxSteps(IdxT& msa) const
    {
      msa = max_steps_;
    }

    template <typename ScalarP, typename IdxP>
    void GridElement<ScalarP, IdxP>::setVariableIndex(IdxT local_index, IdxT global_index)
    {
      variable_indices_[static_cast<size_t>(local_index)] = global_index;
    }

    template <typename ScalarP, typename IdxP>
    IdxP& GridElement<ScalarP, IdxP>::getVariableIndex(IdxT local_index)
    {
      return variable_indices_[static_cast<size_t>(local_index)];
    }

    template <typename ScalarP, typename IdxP>
    void GridElement<ScalarP, IdxP>::setResidualIndex(IdxT local_index, IdxT global_index)
    {
      residual_indices_[static_cast<size_t>(local_index)] = global_index;
    }

    template <typename ScalarP, typename IdxP>
    IdxP& GridElement<ScalarP, IdxP>::getResidualIndex(IdxT local_index)
    {
      return residual_indices_[static_cast<size_t>(local_index)];
    }

    template class GridElement<double, long int>;
    template class GridElement<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit

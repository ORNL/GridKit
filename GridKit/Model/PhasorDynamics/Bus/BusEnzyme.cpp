/**
 * @file BusEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include "BusImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental.
     *
     * This sets values to 0, and these remain unchanged. It is needed to get
     * the indices into the list of entries that will later be deduplicated.
     * Contributions to bus Jacobians from other components are stored in those components.
     *
     * @return int - error code
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Bus..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      if (J_rows_buffer_ == nullptr)
      {
        J_rows_buffer_ = new IdxT[4];
        J_cols_buffer_ = new IdxT[4];
        J_vals_buffer_ = new RealT[4];

        J_rows_buffer_[0] = residual_indices_.at(0);
        J_rows_buffer_[1] = residual_indices_.at(0);
        J_rows_buffer_[2] = residual_indices_.at(1);
        J_rows_buffer_[3] = residual_indices_.at(1);
        J_cols_buffer_[0] = variable_indices_.at(0);
        J_cols_buffer_[1] = variable_indices_.at(1);
        J_cols_buffer_[2] = variable_indices_.at(0);
        J_cols_buffer_[3] = variable_indices_.at(1);
        J_vals_buffer_[0] = 0.0;
        J_vals_buffer_[1] = 0.0;
        J_vals_buffer_[2] = 0.0;
        J_vals_buffer_[3] = 0.0;

        nnz_ = 4;
        this->constructCoo();
      }
      return 0;
    }

    // Available template instantiations
    template class BusBase<double, long int>;
    template class BusBase<double, size_t>;
    template class Bus<double, long int>;
    template class Bus<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit

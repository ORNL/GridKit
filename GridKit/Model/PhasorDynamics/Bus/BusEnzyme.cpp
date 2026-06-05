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
     * @brief Jacobian evaluation experimental. This sets values to 0, for other
     * components to add their contributions.
     *
     * @warning This implementation assumes bus Jacobians are always evaluated
     * _before_ component model Jacobians.
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
        J_.zeroMatrix();

        // Reserve space for a dense matrix of size_*size_.
        // Enyme will compute the appropriate nnz from sparsification.
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
        J_.setValues(1.0, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, 4); //< @todo: Update once sparse storage format changes
      }
      else
      {
        J_.zeroValuedMatrix();
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

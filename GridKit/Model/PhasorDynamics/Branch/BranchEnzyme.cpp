/**
 * @file BranchEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Branch..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      gatherExternalVariables();

      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense blocks.
        // The size of the buffer is the sum of maximum capacities of the blocks.
        // Enyme will compute the appropriate nnz from sparsification.
        auto f_ext_size  = static_cast<size_t>(f_ext_.getSize());
        auto y_ext_size  = static_cast<size_t>(y_ext_.getSize());
        auto buffer_size = f_ext_size * y_ext_size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      GridKit::Enzyme::Sparse::DfDyExt<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::ExternalResidual>::eval(this,
                                                                                                         static_cast<size_t>(f_ext_.getSize()),
                                                                                                         static_cast<size_t>(y_ext_.getSize()),
                                                                                                         residual_indices_ext_.data(),
                                                                                                         variable_indices_ext_.data(),
                                                                                                         y_.getData(),
                                                                                                         yp_.getData(),
                                                                                                         y_ext_.getData(),
                                                                                                         J_rows_buffer_,
                                                                                                         J_cols_buffer_,
                                                                                                         J_vals_buffer_,
                                                                                                         nnz_);

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class Branch<double, long int>;
    template class Branch<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
